/*
    This file reuses a lot of code written by Nicolas Holzschuch
    from "src/subsurface/singlescatter.cpp".

    All code written by me will be marked with [MF2020].
    Unmarked code is by default considered to be written by NH.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sse.h>
#include <mitsuba/core/ssemath.h>
#include "../medium/materials.h"

MTS_NAMESPACE_BEGIN

//////////////////////////////////////////////////////////////////////////////
/// \brief Evaluate the Henyey-Greenstein phase function.
///
Spectrum hg(Float cosTheta, const Spectrum &g) {
	Spectrum temp = Spectrum(1) + g * g + 2 * g * cosTheta;
	return INV_FOURPI * (Spectrum(1) - g * g) / (temp * temp.sqrt());
}

static ref<Mutex> mutex = new Mutex;

class SingleScatterOptimized : public Subsurface {
public:
    SingleScatterOptimized(const Properties &props) : Subsurface(props) {
		/* Single scattering strategy: use fast single scatter? (Jensen) */
		m_fastSingleScatter = props.getBoolean("fastSingleScatter", true);

		/* Single scattering: number of samples along the inside ray ? */
		m_fastSingleScatterSamples = props.getInteger("fssSamples", 2);

		/* Single scattering: use shadow rays? */
		/* This flag only makes sense when debugging, i.e. it should always be
		 * true */
		m_singleScatterShadowRays =
			props.getBoolean("singleScatterShadowRays", true);

		/* Single scattering: compute transmittance? */
		/* This flag only makes sense when debugging, i.e. it should always be
		 * true */
		m_singleScatterTransmittance =
			props.getBoolean("singleScatterTransmittance", true);

		/* Single scattering: number of total internal reflexion? */
		m_singleScatterDepth = props.getInteger("singleScatterDepth", 4);

		/* [MF2020] compute D with or without taking ray differentials into account */
		m_distanceCorrection = props.getBoolean("distanceCorrection", true);

        // Get the material parameters
        lookupMaterial(props, m_sigmaS, m_sigmaA, m_g);
    }

	//-------------------------------------------------------------------------
    virtual void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(BSDF)))
            if (!m_BSDF.get())
                m_BSDF = dynamic_cast<BSDF *>(child);
            else
                Log(EError, "SingleScatterOptimized nodes should have a unique BSDF child.");
        else
            Subsurface::addChild(name, child);
    }

    SingleScatterOptimized(Stream *stream, InstanceManager *manager)
        : Subsurface(stream, manager) {
        m_BSDF = static_cast<BSDF *>(manager->getInstance(stream));
        m_sigmaS = Spectrum(stream);
        m_sigmaA = Spectrum(stream);
        m_g = Spectrum(stream);
        m_eta = stream->readFloat();

        // Additions for single scatter
		m_fastSingleScatter = stream->readBool();
		m_fastSingleScatterSamples = stream->readInt();
		m_singleScatterShadowRays = stream->readBool();
		m_singleScatterTransmittance = stream->readBool();
		m_singleScatterDepth = stream->readInt();
        m_distanceCorrection = stream->readBool();  // [MF2020]
        configure();
    }

    virtual ~SingleScatterOptimized() {}

    void bindUsedResources(ParallelProcess *proc) const {}

    void serialize(Stream *stream, InstanceManager *manager) const {
        Subsurface::serialize(stream, manager);
		manager->serialize(stream, m_BSDF.get());
		m_sigmaS.serialize(stream);
		m_sigmaA.serialize(stream);
		m_g.serialize(stream);
		stream->writeFloat(m_eta);

        // Additions for single scatter
		stream->writeBool(m_fastSingleScatter);
		stream->writeInt(m_fastSingleScatterSamples);
		stream->writeBool(m_singleScatterShadowRays);
		stream->writeBool(m_singleScatterTransmittance);
		stream->writeInt(m_singleScatterDepth);        
        stream->writeBool(m_distanceCorrection);    // [MF2020]
    }

    //-------------- Begin set of functions for single scattering -------------
	Spectrum attenuation(const Spectrum &muT, Float negDistance) const {
		Spectrum result(1.0);
		for (int c = 0; c < SPECTRUM_SAMPLES; ++c)
			if (m_sigmaT[c])
				result[c] = math::fastexp(muT[c] * negDistance);
		return result;
	}

	//------------------------------------------------------------------------
	Spectrum contributionFromThatPoint(const Vector &paramP, const Vector &dPdu,
									   const Vector &dPdv, const Vector &dNsdu,
									   const Vector &dNsdv, const Point &L,
									   const Point &V0, const Vector &dInternal,
									   const Point &P0, const Vector tN[3],
									   const Vector &Ng, Float a11, Float a12,
									   Float a22, const Spectrum &inputSpectrum,
									   const Scene *scene,
									   Float time = 0.0f) const {
		const Point Pc = P0 + paramP[1] * dPdu + paramP[2] * dPdv;
		const Point V = V0 + paramP[0] * dInternal;
		Vector omegaV = V - Pc;
		Float domegaV = omegaV.length();
		omegaV /= domegaV;
		Vector omegaL = L - Pc;
		Float domegaL = omegaL.length();
		omegaL /= domegaL;
		Spectrum result = inputSpectrum;

		if (m_singleScatterShadowRays) {
			// Shadow test 1: is the outgoing point visible from the light
			// source?
			const Ray shadow1 = Ray(Pc, omegaL, ShadowEpsilon,
									domegaL * (1 - ShadowEpsilon), time);
			if (scene->rayIntersect(shadow1)) {
				return Spectrum(0.0f);
			}
		}

		Vector Ns = (tN[0] + paramP[1] * dNsdu + paramP[2] * dNsdv);
		Float idNs = 1.0f / Ns.length();
		Ns *= idNs;
		const Float cosThetaL = dot(omegaL, Ns);
		const Float cosThetaV = dot(omegaV, Ns);

		/* Fresnel transmittance at the new position */
		const Float F = fresnelDielectricExt(cosThetaL, m_eta);

		/* Evaluate the Henyey-Greenstein model */
		const Float cosThetaInternal = dot(omegaV, dInternal);
		Spectrum phase = hg(cosThetaInternal,
							m_g); // reproduces results with +cosThetaInternal.
		result *= (1 - F) * phase;
		result *= m_sigmaS * attenuation(m_sigmaT, -(paramP[0] + domegaV));
		// computing D with ray differentials as in [Walter 2009]
		// u_p and u_s are omega'_V
		const Float mu = cosThetaL + m_eta * cosThetaV;
		// u_p : perpendicular vector
		const Vector u_p = normalize(cross(omegaV, Ns));
		const Vector dPdu_p =
			domegaV * (u_p - (dot(u_p, Ng) / dot(omegaV, Ng)) * omegaV);

		// Normal derivatives are OK
		const Float dudu_p =
			(a22 * dot(dPdu_p, dPdu) - a12 * dot(dPdu_p, dPdv));
		const Float dvdu_p =
			(-a12 * dot(dPdu_p, dPdu) + a11 * dot(dPdu_p, dPdv));
		const Float dwdu_p = -dudu_p - dvdu_p;
		const Vector dNdu_p = dwdu_p * tN[0] + dudu_p * tN[1] + dvdu_p * tN[2];
		const Vector dNndu_p = dNdu_p * idNs - dot(Ns, dNdu_p * idNs) * Ns;
		const Float dmudu_p =
			m_eta * (mu / cosThetaL) * (dot(-u_p, Ns) + dot(omegaV, dNndu_p));
		const Vector domegaLdu_p = m_eta * u_p + dmudu_p * Ns + mu * dNndu_p;
		const Vector L_p =
			dPdu_p - dot(dPdu_p, omegaL) * omegaL + domegaL * domegaLdu_p;

		// u_s : parallel vector
		const Vector u_s = normalize(cross(u_p, omegaV));
		const Vector dPdu_s =
			domegaV * (u_s - (dot(u_s, Ng) / dot(omegaV, Ng)) * omegaV);
		// Normal derivatives
		const Float dudu_s =
			(a22 * dot(dPdu_s, dPdu) - a12 * dot(dPdu_s, dPdv));
		const Float dvdu_s =
			(-a12 * dot(dPdu_s, dPdu) + a11 * dot(dPdu_s, dPdv));
		const Float dwdu_s = -dudu_s - dvdu_s;
		const Vector dNdu_s = dwdu_s * tN[0] + dudu_s * tN[1] + dvdu_s * tN[2];
		const Vector dNndu_s = dNdu_s * idNs - dot(Ns, dNdu_s * idNs) * Ns;
		const Float dmudu_s =
			m_eta * (mu / cosThetaL) * (dot(-u_s, Ns) + dot(omegaV, dNndu_s));
		const Vector domegaLdu_s = m_eta * u_s + dmudu_s * Ns + mu * dNndu_s;
		const Vector L_s =
			dPdu_s - dot(dPdu_s, omegaL) * omegaL + domegaL * domegaLdu_s;
		Float D = cross(L_p, L_s).length();

		// For debug/explanation only: result without the ray-differentials

		if (!m_distanceCorrection) {
		  D = (domegaV + m_eta * domegaL) * (std::abs(cosThetaL/cosThetaV)*domegaV +
											 std::abs(cosThetaV/cosThetaL)*m_eta*domegaL);
		}

		return result / D;
	}

	//-------------------------------------------------------------------------
    Spectrum LoSingle(const Scene *scene, Sampler *sampler,
					  const Intersection &its, const Vector &dInternal,
					  int depth, Float z0) const {
		Spectrum result(0.0f);
		if (depth >= m_singleScatterDepth) {
			return result;
		}

		//---- Look at the intersection
		Ray forwardRay(its.p, dInternal, its.time);
		Intersection its2;
		if (EXPECT_NOT_TAKEN(!scene->rayIntersect(forwardRay, its2))) {
			return result; // starting point seems to be outside the object
		}

		// How large is the object?
		const Float thickness = its2.t;

		// 2014-04-22 jDG: Beware of intersecting mesh: we do have a (possible)
		//                 exit transmittance ONLY when hitting the SAME object
		//                 again.
		if (m_singleScatterTransmittance && its2.shape == getShapes()[0]) {
			// 2014-03-10 jDG: Better to use the <<canonic>> way to get
			//                 both the refracted direction AND the attenuation.
			BSDFSamplingRecord bRec(its2, sampler);
			bRec.typeMask = BSDF::EDeltaTransmission;
			Spectrum bsdfAtt = m_BSDF->sample(bRec, sampler->next2D());
			sampler->advance();
			Vector dOutgoing = its2.toWorld(bRec.wo);

			// compute radiance through the object
			if (!bsdfAtt.isZero()) {
				// There is transmitted light ?
				// 2014-04-22 jDG: added indirect surface radiance (so take ALL
				// of them)
				RadianceQueryRecord query(scene, sampler);
				query.newQuery(
					RadianceQueryRecord::ERadiance // All types of radiances
												   // (including us) !
						| RadianceQueryRecord::EIntersection, // Make sure to
															  // compute
															  // intersect
															  // first.
					its.shape->getExteriorMedium());
				query.depth = depth;
				Spectrum refracted = m_integrator->Li(
					RayDifferential(its2.p, dOutgoing, its.time), query);
				refracted *= bsdfAtt;
				refracted *= attenuation(m_sigmaT, -thickness);
				result += refracted;
			}
		}

		// Internal reflection (total or quasi total). Reccursivly call LoSingle
		// again.
		{
			BSDFSamplingRecord bRec(its2, sampler);
			const BSDF *bsdf = 0;
			if (its2.shape == getShapes()[0]) {
				bRec.typeMask =
					BSDF::EReflection; // Evrything except exit rays.
				bsdf = m_BSDF.get();
			} else {
				bRec.typeMask = BSDF::EAll; // Everything.
				bsdf = its2.getBSDF();
			}

			Spectrum bsdfAtt = bsdf->sample(bRec, sampler->next2D());
			sampler->advance();
			Vector dReflect = its2.toWorld(bRec.wo);

			if (!bsdfAtt.isZero()) {
				Spectrum reflected = LoSingle(scene, sampler, its2, dReflect,
											  depth + 1, z0 + thickness);
				reflected *= bsdfAtt;
				// 2014-02-24 jDG: Attenuation on the whole path, at once.
				reflected *= attenuation(m_sigmaT, -thickness);
				result += reflected;
			}
		}

		/* Sample a point on a light source */
		DirectSamplingRecord dRec(its.p, its.time);
		// 2014-04-30 NH: Added the m_eta^2 coefficient, for the light ray enter
		// the material.
		const Spectrum value =
			m_eta * m_eta *
			scene->sampleEmitterDirect(dRec, sampler->next2D(), false);
		sampler->advance();
		if (value.isZero())
			return result;

		const Point L = dRec.p;
		if (m_fastSingleScatter) {
			// Classical SS approximation. Shoot one ray back to the light
			// source.
			// We allow more than one sample along the ray, though.
			const Float sMax = 1 - exp(-(thickness / m_radius));
			const Float dSamples =
				m_fastSingleScatterSamples
					? sMax / Float(m_fastSingleScatterSamples)
					: sMax;

			const Spectrum weight0 =
				(dSamples * m_radius * (dRec.dist * dRec.dist)) * m_sigmaS;

			for (int s = 0; s < m_fastSingleScatterSamples; ++s) {
				Float sample = sampler->next1D() * sMax;
				sampler->advance();
				const Float dist = -math::fastlog(1 - sample) * m_radius;
				const Point V = its.p + dist * dInternal;
				if (EXPECT_NOT_TAKEN(dist > thickness))
					continue;

				/* Sampling weight, because sampling should have been
				 * chanel-dependant */
				Spectrum weight = weight0 * math::fastexp(m_invRadius * dist);

				/* First, connect to the light source */
				Vector VL = L - V;
				Float dVL = VL.length();
				VL /= dVL;
				Ray toTheLight(V, VL, Epsilon, dVL * (1 - ShadowEpsilon), its.time);
				if (!scene->rayIntersect(toTheLight, its2))
					continue;

				const Point PWorld = its2.p;
				/* Make sure that the light source is not occluded from this
				 * position */
				Vector omegaL = L - PWorld;
				Float dL = omegaL.length();
				omegaL /= dL;

				/* shadow ray */
				Ray ray(PWorld, omegaL, Epsilon, dL * (1 - ShadowEpsilon),
						its2.time);
				if (scene->rayIntersect(ray)) {
					continue;
				}

				Vector omegaV = V - PWorld;
				Float dV = omegaV.length();
				omegaV /= dV;

				/* Account for importance sampling wrt. transmittance */
				const Float cosThetaL = dot(omegaL, its2.shFrame.n);
				const Float cosThetaV = dot(omegaV, its2.shFrame.n);
				if (cosThetaL == 0 || cosThetaV == 0)
					continue;

				/* Fresnel transmittance at the new position */
				const Float F = fresnelDielectricExt(cosThetaL, m_eta);

				/* Evaluate the Henyey-Greenstein model */
				Float cosThetaInternal = dot(omegaV, dInternal);
				Spectrum phase;
				phase = hg(cosThetaInternal, m_g) *
						attenuation(m_sigmaT, -(dist + dV));

				const Float D = (dV + m_eta * dL) *
								(std::abs(cosThetaL / cosThetaV) * dV +
								 std::abs(cosThetaV / cosThetaL) * m_eta * dL);
				result += ((1 - F) / D) * phase * value * weight;
			}
		}
    }

    //---------------- End set of functions for single scattering ----------------------
    Spectrum Lo(const Scene *scene, Sampler *sampler, const Intersection &its,
                const Vector &d, int depth) const {
        //---- Initialize intergator and BSDF stuff from the first intersection seen.
		if (!m_integrator) {
			LockGuard lock(mutex);
			if (!m_integrator) {
				SingleScatterOptimized *self = const_cast<SingleScatterOptimized *>(this);

				self->m_integrator = dynamic_cast<const MonteCarloIntegrator *>(
					scene->getIntegrator());
				if (!m_integrator)
					Log(EError, "Single scatter requires a sampling-based "
								"surface integrator!");
				if (!m_integrator->getClass()->derivesFrom(
						MTS_CLASS(SamplingIntegrator)))
					Log(EError, "Single scatter requires a sampling-based "
								"surface integrator!");
			}
		}

		Spectrum result(0.0f);

		//---- Perform Reflections (if any) ----------------------------------
		{
			BSDFSamplingRecord bRec(its, sampler);
			bRec.typeMask = BSDF::EDeltaReflection;
			Spectrum reflectAttenuation =
				m_BSDF->sample(bRec, sampler->next2D());
			Vector dBounced = its.toWorld(bRec.wo);
			sampler->advance();
			if (!reflectAttenuation.isZero()) {
				RadianceQueryRecord query(scene, sampler);
				query.newQuery(
					RadianceQueryRecord::ERadiance // All radiances sources.
						| RadianceQueryRecord::EIntersection, // Should compute
															  // first
															  // intersect.
					its.shape->getExteriorMedium());
				query.depth = depth + 1;
				RayDifferential ray(its.p, dBounced, its.time);
				result += reflectAttenuation * m_integrator->Li(ray, query);
			}
		}
		//---- Perform refractions (if any) and single scatter ----------------------
		{
			BSDFSamplingRecord bRec(its, sampler);
			bRec.typeMask = BSDF::EDeltaTransmission;
			Spectrum refractAttenuation =
				m_BSDF->sample(bRec, sampler->next2D());
			Vector dInternal = its.toWorld(bRec.wo);
			sampler->advance();

			if (!refractAttenuation.isZero()) {
				result +=
					refractAttenuation *
					LoSingle(scene, sampler, its, dInternal, depth + 1, 0);
			}
		}
		return result;   
    }

    void configure() {
		if (!m_BSDF.get()) {
			Log(EError, "Single scatter should have a BSDF child node.");
			m_eta = 1;
		} else
			m_eta = m_BSDF->getEta();
		m_invEta = 1 / m_eta;

		if (m_eta < 1)
			Log(EError, "Unsupported material configuration (intIOR/extIOR < 1)");

		m_sigmaT = m_sigmaA + m_sigmaS;

		/* Find the smallest mean-free path over all wavelengths */
		Spectrum mfp = Spectrum(1.0f) / m_sigmaT;
		m_radius = std::numeric_limits<Float>::max();
		for (int lambda = 0; lambda < SPECTRUM_SAMPLES; lambda++)
			m_radius = std::min(m_radius, mfp[lambda]);
		m_invRadius = 1.0f / m_radius;

    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
                    const RenderJob *job, int sceneResID, int cameraResID,
                    int _samplerResID) {
		if (!scene->getIntegrator()->getClass()->derivesFrom(
				MTS_CLASS(SamplingIntegrator)))
			Log(EError, "The single scattering pluging requires "
						"a sampling-based surface integrator!");
		return true;

    }

    void wakeup(ConfigurableObject *parent,
                std::map<std::string, SerializableObject *> &params) {}

    void cancel() {}

    MTS_DECLARE_CLASS()
private:
    ref<const MonteCarloIntegrator> m_integrator;

    Float m_radius, m_invRadius;    // smallest mean-free path over all wavelengths
    Float m_eta, m_invEta;          // relative index of refraction (intIOR/extIOR)
    Spectrum m_sigmaS, m_sigmaA, m_sigmaT, m_g; // medium properties
    ref<const BSDF> m_BSDF;

    bool m_fastSingleScatter;
	int m_fastSingleScatterSamples;
	bool m_singleScatterShadowRays;
	bool m_singleScatterTransmittance;
	int m_singleScatterDepth;
    bool m_distanceCorrection;  // [MF2020]
};

MTS_IMPLEMENT_CLASS_S(SingleScatterOptimized, false, Subsurface)
MTS_EXPORT_PLUGIN(SingleScatterOptimized, "Single scattering model (new version)");
MTS_NAMESPACE_END