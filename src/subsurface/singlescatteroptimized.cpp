/*
    This file reuses a lot of code written by Nicolas Holzschuch
    from "src/subsurface/singlescatter.cpp".

	Look for flags ([TODO], [CHECK])
	
	[TODO] change relevant failcases to asserts
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sse.h>
#include <mitsuba/core/ssemath.h>
#include "../medium/materials.h"

MTS_NAMESPACE_BEGIN

//-------------- [MF2020] Begin structures ------------------------------------

enum EConicType {
	EEllipse = 0,
	EHyperbola,
	EParabola,
	EDegenerate
};

enum EConicBranch {
	EPositive = 0,
	ENegative,
	ESingle		// Ellipses have a single branch
};

enum ETriangleSide {
	ESideA = 0,	// Horizontal side, b = 0
	ESideB,		// Vertical side, a = 0
	ESideAB		// Diagonal side, a + b = 1
};

struct ConicParam {
	Float value;			// parametrisation in canonical conic space
	EConicBranch branch;	// +/- branch for hyperbolas
};

// Contains the barycentric and parametric representation of a point on the conic section
struct ConicPoint {
	Vector3 position;		// homogeneous coordinates in barycentric space (a, b, 1)
	ConicParam parameter;
};

// Intersection points between the conic section and the triangle
struct IsectPoint {
	ConicPoint cPoint;
	ETriangleSide side;		// edge of the triangle the intersection point lies on
};

struct ConicChord {
	ConicPoint q0;		// start point of the chord
	ConicPoint q1;		// end point of the chord
	ConicPoint qMid;	// conic point maximising orthogonal distance to the chord
	Float distToConic;	// distance to the chord
	std::vector<Float> polynomial;
};

struct ConicSection {
	// Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
	Float coefficients[6];	// A through F correspond to 0 through 5
	Float k;				// = detConic / detReduced
	EConicType conicType;
	Vector2 center;
	Float semiMajorAxis;	// positive
	Float semiMinorAxis;	// positive
	Float angle;			// in [-pi, pi]
	Matrix3x3 transformCanonicalToBarycentric;	// transform matrices
	Matrix3x3 transformBarycentricToCanonical;
};

//-------------- [MF2020] End structures ------------------------------------

//-------------- [MF2020] Begin utility functions -----------------------------

/*	Comparison function to sort intersection points according to parameter value and branch
	
	Lexicographic order on pairs consisting of (branch, value), where the order on 
	values is the classical order on real values and the ordering on branches is
	brNeg < brPos < brSingle
	[Declared outside of class because of compiler error: 'invalid use of non-static member function']
*/
bool compIsectPts(IsectPoint iPt1, IsectPoint iPt2) {
	// Branch order
	if ((iPt1.cPoint.parameter.branch == EConicBranch::ENegative) &&
		(iPt2.cPoint.parameter.branch != EConicBranch::ENegative)) {
		return true;		
	}
	if ((iPt1.cPoint.parameter.branch == EConicBranch::EPositive) &&
		(iPt2.cPoint.parameter.branch == EConicBranch::ESingle)) {
		return true;
	}

	// Value order
	return iPt1.cPoint.parameter.value < iPt2.cPoint.parameter.value;
}

// Taken from Mitsuba code, tests for equality with zero
bool isAlmostZero(Float x) {
	return std::abs(x) < RCPOVERFLOW;
}

int sign(Float val) {
	return (Float(0) < val) - (val < Float(0));
}

bool haveDiffSigns(Float a, Float b) {
	return sign(a) * sign(b) == -1;
}

int signChanges(const std::vector<Float>& seq) {
	int size = seq.size();
	if (size == 0) return 0;

	int firstNonZero = 0;
	while (isAlmostZero(seq[firstNonZero])) ++firstNonZero;

	int numChanges = 0;
	int i = firstNonZero + 1;
	int lastNonZero = firstNonZero;
	while (i < size) {
		if (!isAlmostZero(seq[i])) {
			if (sign(seq[lastNonZero] * seq[i]) == -1) ++numChanges;
			lastNonZero = i;
		}
		++i;
	}

	return numChanges;
}

//-------------- [MF2020] End utility functions -----------------------------

//-------------- [MF2020] Begin set of functions for polynomials -----------

/*	poly[i] represents the polynomial corresponding to the monomial X^i
	The degree of the polynomial can be obtained by substracting 1 to the size
	of the underlying vector
	[TODO] Replace with Boost polynomials?
	[CHECK] The assumption is made that degree == size - 1
		(does not account for null leading coefficients)
*/

void polynomialPrint(const std::vector<Float>& poly) {
	for (int i = 0; i < poly.size() - 1; ++i) {
		std::cout << poly[i] << "\t";
	}
	std::cout << poly[poly.size() - 1] << std::endl;
}

std::vector<Float> polynomialDerivative(const std::vector<Float>& poly) {
	std::vector<Float> polyDerivative(poly.size() - 1);

	for (int i = 0; i < polyDerivative.size(); ++i) {
		polyDerivative[i] = (i+1) * poly[i+1];
	}

	return polyDerivative;
}

//	Horner evaluation of a polynomial, taken from Higham 2002 (p94)
Float polynomialHorner(const std::vector<Float>& poly, Float t) {
	int n = poly.size() - 1;
	Float result = poly[n];

	for (int i = n - 1; i >= 0; --i) {
		result = t * result + poly[i];
	}

	return result;
}

/*	Horner evaluation of a polynomial and its first derivative, from Higham 2002 (p97, k = 1)
	valP contains the result of evaluating the polynomial,
	valD contains the result of evaluating its first derivative.
*/
void polynomialHornerDeriv(const std::vector<Float>& poly, Float t,
						   Float& valP, Float& valD) {
	int n = poly.size() - 1;
	valP = poly[n];
	valD = 0;

	for (int j = n - 1; j >= 0; --j) {
		valD = t * valD + valP;
		valP = t * valP + poly[j];
	}
}

void polynomialAdd(std::vector<Float>& polyDest, const std::vector<Float>& polySrc) {
	int newSize = std::max(polyDest.size(), polySrc.size());
	polyDest.resize(newSize, 0.0);

	for (int i = 0; i < newSize; ++i) {
		polyDest[i] += polySrc[i];
	}
}

void polynomialMult(std::vector<Float>& polyDest, const std::vector<Float>& polySrc) {
	int degDest = polyDest.size() - 1;
	int degSrc = polySrc.size() - 1;
	int newDeg = degDest + degSrc;
	// newSize = newDeg + 1 and newDeg = degDest + degSrc = sizeDest + sizeSrc - 2
	std::vector<Float> polyRes(newDeg + 1, 0.0);

	for (int k = 0; k <= newDeg; ++k) {
		for (int i = std::max(0, k - degSrc); i <= std::min(k, degDest); ++i){
			polyRes[k] += polyDest[i] * polySrc[k - i];
		}
	}
	polyDest = polyRes;
}

void polynomialSMult(std::vector<Float>& polyDest, Float k) {
	for (int i = 0; i < polyDest.size(); ++i) {
		polyDest[i] *= k;
	}
}

void polynomialNeg(std::vector<Float>& polyDest) {
	for (int i = 0; i < polyDest.size(); ++i) {
		polyDest[i] = -polyDest[i];
	}
}

// Computes the remainder of the division of polyDividend by polyDivisor
std::vector<Float> polynomialRemainder(const std::vector<Float>& polyDividend,
									   const std::vector<Float>& polyDivisor) {
	std::vector<Float> polyRemainder(polyDividend);
	int sizeR = polyRemainder.size();
	int sizeD = polyDivisor.size();
	
	// std::cout << "#################################################" << std::endl;
	// std::cout << "polyDividend: (" << polyDividend.size() << ")\t";
	// polynomialPrint(polyDividend);
	
	// std::cout << "polyDivisor: (" << polyDivisor.size() << ")\t";
	// polynomialPrint(polyDivisor);	

	Float leadingCoeff = polyDivisor[sizeD - 1];
	// std::cout << "leadingCoeff: " << leadingCoeff << std::endl;
	while (sizeR >= sizeD) {
		Float factor = polyRemainder[sizeR - 1] / leadingCoeff;
		// std::cout << "factor: " << factor << std::endl;
		int degDiff = sizeR - sizeD;
		// std::cout << "degDiff: " << degDiff << std::endl;
		
		std::vector<Float> tmpDivisor(sizeR, 0);
		for (int i = 0; i < sizeD; ++i) {
		    tmpDivisor[degDiff+i] = polyDivisor[i];
		}

		// std::cout << "tmpDivisor: (" << tmpDivisor.size() << ")\t";
		// polynomialPrint(tmpDivisor);

		polynomialSMult(tmpDivisor, -factor);

		// std::cout << "tmpDivisor: (" << tmpDivisor.size() << ")\t";
		// polynomialPrint(tmpDivisor);

		polynomialAdd(polyRemainder, tmpDivisor);

		// std::cout << "polyRemainder: (" << polyRemainder.size() << ")\t";
		// polynomialPrint(polyRemainder);

		int leadingZeros = 1;
		for (int i = sizeR - 2; i >= 0; --i) {
			if (!isAlmostZero(polyRemainder[i])) break;
			++leadingZeros;
		}
		sizeR -= leadingZeros;
		// std::cout << "leadingZeros: " << leadingZeros << std::endl;
		// std::cout << "sizeR: " << sizeR << std::endl;
		polyRemainder.resize(sizeR);

		// std::cout << "polyResized: (" << polyRemainder.size() << ")\t";
		// polynomialPrint(polyRemainder);
	}
	
	return polyRemainder;
}

int sturmAlgorithm(const std::vector<Float>& poly, Float a, Float b) {
	Float polyA, polyDerivA, polyB, polyDerivB;
	polynomialHornerDeriv(poly, a, polyA, polyDerivA);
	polynomialHornerDeriv(poly, b, polyB, polyDerivB);

	Float lastNonZeroA = 0;
	Float lastNonZeroB = 0;

	bool hasLNZA = true;
	bool hasLNZB = true;

	int numChangesA = 0;
	int numChangesB = 0;

	/*	Initialisation based on whether evaluating P0 and P1 at a value t
		provides zero or nonzero results
	*/
	if (!isAlmostZero(polyDerivA)) {
		lastNonZeroA = polyDerivA;

		if (!isAlmostZero(polyA)) {
			if (haveDiffSigns(polyA, polyDerivA)) ++numChangesA;
		}
	} else {
		if (!isAlmostZero(polyA)) {
			lastNonZeroA = polyA;
		} else {
			hasLNZA = false;
		}
	}

	if (!isAlmostZero(polyDerivB)) {
		lastNonZeroB = polyDerivB;

		if (!isAlmostZero(polyB)) {
			if (haveDiffSigns(polyB, polyDerivB)) ++numChangesB;
		}
	} else {
		if (!isAlmostZero(polyB)) {
			lastNonZeroB = polyB;
		} else {
			hasLNZB = false;
		}
	}

	std::vector<Float> polyPrev(poly);
	std::vector<Float> polyCurr = polynomialDerivative(poly);
	
	// Algorithm
	while (polyCurr.size() > 1) {
		// Next polynomial in the sequence
		std::vector<Float> polyTmp = polynomialRemainder(polyPrev, polyCurr);
		polyPrev = polyCurr;
		polyCurr = polyTmp;
		polynomialNeg(polyCurr);
		
		Float evalA = polynomialHorner(polyCurr, a);
		if (!hasLNZA) {		// No nonzero result yet
			if (!isAlmostZero(evalA)) {
				hasLNZA = true;
				lastNonZeroA = evalA;
			}
		} else {
			// Compare sign of current eval (if nonzero) to last nonzero eval
			if (!isAlmostZero(evalA)) {
				if (haveDiffSigns(lastNonZeroA, evalA)) {
					++numChangesA;
					lastNonZeroA = evalA;
				}
			}
		}

		Float evalB = polynomialHorner(polyCurr, b);
		if (!hasLNZB) {		// No nonzero result yet
			if (!isAlmostZero(evalB)) {
				hasLNZB = true;
				lastNonZeroB = evalB;
			}
		} else {
			// Compare sign of current eval (if nonzero) to last nonzero eval
			if (!isAlmostZero(evalB)) {
				if (haveDiffSigns(lastNonZeroB, evalB)) {
					++numChangesB;
					lastNonZeroB = evalB;
				}
			}
		}
	}

	return numChangesA - numChangesB;
}

/*	Generates the quadratic term in the constraint function
	\norm{X - P(t)}^2
	where the variable t is the position along the chord
	poly must be of size 3 (deg 2)
*/
void polynomialQuadraticConstraint(std::vector<Float>& poly, Vector3 X,
								   Vector3 Pc, Vector3 Pt) {
	Vector3 XPc = X - Pc;
	poly[0] = XPc.lengthSquared();
	poly[1] = -2.0 * dot(XPc, Pt);
	poly[2] = Pt.lengthSquared();
}

/*	Generates the quartic term in the constraint function
	\norm{N(t) \wedge (X - P(t))}^2
	where the variable t is the position along the chord
	poly must be of size 5 (deg 4)
*/
void polynomialQuarticConstraint(std::vector<Float>& poly, Vector3 X,
								 Vector3 Pc, Vector3 Pt,
								 Vector3 Nc, Vector3 Nt) {
	Vector3 XPc = X - Pc;
	Vector3 A = cross(Nc, XPc);
	Vector3 B = cross(Nt, XPc) - cross(Nc, Pt);
	Vector3 C = cross(Nt, Pt);

	poly[0] = C.lengthSquared();
	poly[1] = 2.0 * dot(B, C);
	poly[2] = B.lengthSquared() + 2.0 * dot(A, C);
	poly[3] = 2.0 * dot(A, B);
	poly[4] = A.lengthSquared();
}

int polynomialDescartes(const std::vector<Float>& poly) {
	std::vector<Float> polyTransformed(poly);

	int n = poly.size() - 1;
	for (int i = 1; i <= n; ++i) {
		for (int k = i; k >= 1; --k) {
			polyTransformed[k] += polyTransformed[k+1];
		}
	}

	return signChanges(polyTransformed);
}

/* 	Naïve Newton's method
	[TODO] Improve method robustness
*/
bool polynomialNewton(const std::vector<Float>& poly, Float& root) {
	int maxIter = 10;
	Float minPrecision = 1e-10;
	int iter = 0;
	bool foundSol = false;

	Float val0 = polynomialHorner(poly, 0.0);
	Float val1 = polynomialHorner(poly, 1.0);
	// root of the linear interpolation of P in [0, 1]
	Float x = val0 / (val0 - val1);
	x = math::clamp(x, (Float)0.0, (Float)1.0);	// [CHECK] is this necessary?

	Float valP, valD;
	polynomialHornerDeriv(poly, x, valP, valD);
	if (std::abs(valP) < minPrecision) {	// check if initial guess is a solution
		foundSol = true;
	}

	// Newton's method, stops after a number of iterations
	// or when the step falls below a threshold
	while (!foundSol && (iter < maxIter)) {
		Float step = - valP / valD;
		Float oldX = x;
		x += step;
		++iter;

		if (std::abs(x - oldX) < minPrecision) {
			foundSol = true;
		}

		polynomialHornerDeriv(poly, x, valP, valD);
		// std::cout << "iter: " << iter << "\tx = " << x << "\tvalP = " << valP << std::endl;
	}

	if (foundSol) {	// valid solution only if between 0 and 1
		if ((x >= 0.0) && (x <= 1.0)) root = x;
		else foundSol = false;
	}
	
	return foundSol;
}

//-------------- [MF2020] End set of functions for polynomials -------------

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

		/* [MF2020] Single scattering: number of samples along the internal ray */
		m_singleScatterSamples = props.getInteger("singleScatterSamples", 1);

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
		m_singleScatterSamples = stream->readInt();	// [MF2020]
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
		stream->writeInt(m_singleScatterSamples);	// [MF2020]
    }

	//-------------- [MF2020] Begin set of functions for conic sections -------
	bool isInTriangle(Float a, Float b) const {
		return (a >= 0.0 && a <= 1.0) &&
			   (b >= 0.0 && b <= 1.0) &&
			   (a + b <= 1.0);
	}

	bool computeConicParameters(ConicSection &conic) const {
		Matrix3x3 conicMatrix(conic.coefficients[0], 0.5f * conic.coefficients[1], 0.5f * conic.coefficients[3],
							  0.5f * conic.coefficients[1], conic.coefficients[2], 0.5f * conic.coefficients[4],
							  0.5f * conic.coefficients[3], 0.5f * conic.coefficients[4], conic.coefficients[5]);
		Matrix2x2 conicReduced(conic.coefficients[0], 0.5f * conic.coefficients[1],
							   0.5f * conic.coefficients[1], conic.coefficients[2]);

		Float detConic = conicMatrix.det();
		Float detReduced = conicReduced.det();
		conic.k = detConic / detReduced;

		// Compute the type of the coplanarity conic
		// [TODO] Use isAlmostZero
		if(std::abs(detConic) < RCPOVERFLOW) {
			conic.conicType = EConicType::EDegenerate;
		} else {
			if (std::abs(detReduced) <= RCPOVERFLOW) {
				conic.conicType = EConicType::EParabola;
			} else if (detReduced > 0.0) {
				conic.conicType = EConicType::EEllipse;
			} else {
				conic.conicType = EConicType::EHyperbola;
			}
		}

		// Fails if the conic is not a central conic, or is degenerate
		if ((conic.conicType == EConicType::EParabola) ||
			(conic.conicType == EConicType::EDegenerate)) {
			if (conic.conicType == EConicType::EDegenerate)
				Log(EError, "computeConicParameters() failed: 1a");
			else
				Log(EError, "computeConicParameters() failed: 1b");
			
			return false;
		}

		// Compute the CENTER, the ANGLE and the SEMI-MINOR/MAJOR AXES of the central conic
		// Only works for central conics (Hyperbolas/Ellipses)
		// Compute the center of the conic
		const Float a[2][2] = {{conicReduced(0, 0), conicReduced(0, 1)},
							   {conicReduced(1, 0), conicReduced(1, 1)}};
		const Float b[2] = {-0.5f * conic.coefficients[3],
							-0.5f * conic.coefficients[4]};
		Float x[2] = {0};
		
		// Passes if conic is not a parabola
		if (!solveLinearSystem2x2(a, b, x))
			return false;
		
		conic.center[0] = x[0];
		conic.center[1] = x[1];

		// Compute the angle of the conic
		Float eigval0, eigval1;	// Reduced matrix eigenvalues
		// Should never fail
		if (!solveQuadratic(1.0f, -(conic.coefficients[0] + conic.coefficients[2]),
							detReduced, eigval0, eigval1)) {
			Log(EError, "computeConicParameters() failed: 2");
			return false;
		}

		// Compute the semi-minor/major axes of the conic
		Float majorEigval = 0.0f; 
		Float minorEigval = 0.0f;
		bool condEigval;

		switch (conic.conicType) {
			case EConicType::EEllipse:
				// the eigenvalue with the smaller magnitude corresponds to the major axis
				condEigval = std::abs(eigval0) < std::abs(eigval1);
				majorEigval = condEigval ? eigval0 : eigval1;
				minorEigval = condEigval ? eigval1 : eigval0;
				break;
			case EConicType::EHyperbola:
				// the eigenvalue with a sign opposite to conic.k corresponds to the transverse axis
				condEigval = eigval0 * conic.k < 0;
				majorEigval = condEigval ? eigval0 : eigval1;
				minorEigval = condEigval ? eigval1 : eigval0;
				break;
			default:
				break;
		}

		conic.semiMajorAxis = std::pow(std::abs(conic.k / majorEigval), 0.5);
		conic.semiMinorAxis = std::pow(std::abs(conic.k / minorEigval), 0.5);

		conic.angle = std::atan(2.0f * (majorEigval - conic.coefficients[0]) / conic.coefficients[1]);

		// Compute transformation matrix between barycentric and canonical conic space
		// (canonical --> where the conic is either a circle or a unit hyperbola)
		const Matrix3x3 translation(1.0, 0.0, conic.center[0],
							  0.0, 1.0, conic.center[1],
							  0.0, 0.0, 1.0);
		const Matrix3x3 rotation(std::cos(conic.angle), -std::sin(conic.angle), 0.0,
						   std::sin(conic.angle), std::cos(conic.angle), 0.0,
						   0.0, 0.0, 1.0);
		const Matrix3x3 scaling(conic.semiMajorAxis, 0.0, 0.0,
						  0.0, conic.semiMinorAxis, 0.0,
						  0.0, 0.0, 1.0);
		
		Matrix3x3 transform = translation * rotation * scaling;
		conic.transformCanonicalToBarycentric = transform;
		Matrix3x3 transformInv(0.0);
		if (!transform.invert(transformInv)){
			Log(EError, "computeConicParameters() failed: 3");
			return false;
		}
		conic.transformBarycentricToCanonical = transformInv;

		return true;
	}

	/*	Takes the barycentric coordinates of a point on a conic section and returns
		the corresponding parameter in the conic canonical space.
	*/
	ConicParam conicPosToParam(const ConicSection &conic, const Vector position) const {
		Vector3 canonicalPos = conic.transformBarycentricToCanonical * position;
		
		ConicParam cParam = {0.0, EConicBranch::ESingle};
		if (conic.conicType == EConicType::EEllipse) {
			cParam.value = std::atan2(canonicalPos[1], canonicalPos[0]);
		} else if (conic.conicType == EConicType::EHyperbola) {
			cParam.value = std::asinh(canonicalPos[1]);
			cParam.branch = (canonicalPos[0] > 0) ? EConicBranch::EPositive :EConicBranch::ENegative;
		}

		return cParam;
	}

	/*	Takes the parameter corresponding to a point in the conic canonical space
		and returns its corresponding barycentric coordinates.
	*/
	Vector3 conicParamToPos(const ConicSection &conic, const ConicParam cParam) const {
		Vector3 canonicalPos(0.0);
		canonicalPos[2] = 1.0;	// Homogeneous coordinates

		if (conic.conicType == EConicType::EEllipse) {
			canonicalPos[0] = std::cos(cParam.value);
			canonicalPos[1] = std::sin(cParam.value);
		} else if (conic.conicType == EConicType::EHyperbola) {
			canonicalPos[0] = std::cosh(cParam.value);
			if (cParam.branch == EConicBranch::ENegative)
				canonicalPos[0] *= -1.0;
			canonicalPos[1] = std::sinh(cParam.value);
		}

		Vector3 barycentricPos = conic.transformCanonicalToBarycentric * canonicalPos;

		return barycentricPos;
	}

	/*	Computes the best point on the conic for subdividing the current chord
		and the maximal distance from the chord to that point.
	*/
	bool chordSplit(const ConicSection& conic, ConicChord& chord) const {
		Float qMidParamValue = 0.5 * (chord.q0.parameter.value + chord.q1.parameter.value);
		ConicParam qMidParam = { .value = qMidParamValue, .branch = chord.q0.parameter.branch };
		Vector3 qMidPosition = conicParamToPos(conic, qMidParam);
		
		if (!isInTriangle(qMidPosition[0], qMidPosition[1])) {
			// Should only happen for ellipses
			if (conic.conicType != EConicType::EEllipse) return false;
			qMidParam.value += M_PI;
			qMidPosition = conicParamToPos(conic, qMidParam);
			// One of the two points should always be inside the triangle
			if (!isInTriangle(qMidPosition[0], qMidPosition[1])) return false;
		}
	
		chord.qMid.position = qMidPosition;
		chord.qMid.parameter = qMidParam;

		Vector3 q1q0 = chord.q1.position - chord.q0.position;
		Vector3 qMidq0 = chord.qMid.position - chord.q0.position;

		chord.distToConic = cross(q1q0, qMidq0).length() / q1q0.length();

		return true;
	}

	/*	Generates the whole constraint polynomial over a given chord
	If we denote the quadratic terms by Q_2(X) and the quartic terms by Q_4(X),
	the constraint polynomial is defined by:
	eta^2 * Q_2(L) * Q_4(V) - Q_2(V) * Q_4(L)
	*/
	void chordPolynomial(ConicChord& chord, Float eta,
					     Vector P10, Vector P20, Vector vP0,
						 Vector N10, Vector N20, Vector vN0,
						 Vector vL, Vector vV) const {
		Vector3 Delta = chord.q1.position - chord.q0.position;

		Matrix3x3 pBar(P10, P20, vP0);
		Vector3 Pc = pBar * chord.q0.position;
		Vector3 Pt = pBar * Delta;

		Matrix3x3 nBar(N10, N20, vN0);
		Vector3 Nc = nBar * chord.q0.position;
		Vector3 Nt = nBar * Delta;

		// LHS stored in chord.polynomial
		chord.polynomial.resize(3);	// std::vector<Float> polyQuadraticLHS(3);
		std::vector<Float> polyQuarticLHS(5);
		polynomialQuadraticConstraint(chord.polynomial, vL, Pc, Pt);
		polynomialQuarticConstraint(polyQuarticLHS, vV, Pc, Pt, Nc, Nt);
		polynomialMult(chord.polynomial, polyQuarticLHS);	// resize to 7
		polynomialSMult(chord.polynomial, eta * eta);

		// RHS stored in polyQuadraticRHS
		std::vector<Float> polyQuadraticRHS(3);
		std::vector<Float> polyQuarticRHS(5);
		polynomialQuadraticConstraint(polyQuadraticRHS, vV, Pc, Pt);
		polynomialQuarticConstraint(polyQuarticRHS, vV, Pc, Pt, Nc, Nt);
		polynomialMult(polyQuadraticRHS, polyQuarticRHS);	// resize to 7
		polynomialSMult(polyQuadraticRHS, -1.0);

		// Full constraint polynomial stored in chord.polynomial
		polynomialAdd(chord.polynomial, polyQuadraticRHS);
	}

	/*	Creates the subdivision of the conic arc by recursively splitting each chord in two
		until every chord is close enough to the conic arc
	*/
	// [CHECK]
	// chords' elements should already have midpoint information (call chordSplit beforehand)
	void conicApproximationSubdivide(std::list<ConicChord>& chords,
									 const ConicSection& conic, Float minPrecision) const {
		std::list<ConicChord>::iterator it = chords.begin();
		
		while (it != chords.end()) {
			if (it->distToConic > minPrecision) {
				ConicChord chordLeft = { .q0 = it->q0, .q1 = it->qMid };
				ConicChord chordRight = { .q0 = it->qMid, .q1 = it->q1 };
				chordSplit(conic, chordLeft);
				chordSplit(conic, chordRight);
				chords.push_back(chordLeft);
				chords.push_back(chordRight);

				it = chords.erase(it);	// Returns iterator pointing to next element
			} else {
				++it;
			}
		}
	}

	/*	Discards the conic chords whose polynomial doesn't have any roots in [0, 1]
		Uses Descartes' rule for estimating the number of roots of a polynomial
		[TODO] Implement true VCA algorithm
	*/
	void conicApproximationDescartes(std::list<ConicChord>& chords, const ConicSection& conic,
								 Vector P10, Vector P20, Vector vP0,
						 		 Vector N10, Vector N20, Vector vN0,
								 Vector vL, Vector vV) const {
		std::list<ConicChord>::iterator it = chords.begin();
		
		while (it != chords.end()) {
			chordPolynomial(*it, m_eta, P10, P20, vP0, N10, N20, vN0, vL, vV);
			int descartesBound = polynomialDescartes(it->polynomial);
			//std::cout << "descartesBound = " << descartesBound << std::endl;
			if (descartesBound == 0) {	// No roots, discard the chord
				it = chords.erase(it);
			} else if (descartesBound == 1) { // Single root, keep the chord
				++it;
			} else {	// More than one root over the chord
				// ConicChord chordLeft = { .q0 = it->q0, .q1 = it->qMid };
				// ConicChord chordRight = { .q0 = it->qMid, .q1 = it->q1 };
				// chordSplit(conic, chordLeft);
				// chordSplit(conic, chordRight);
				// chords.push_back(chordLeft);
				// chords.push_back(chordRight);
				it = chords.erase(it);	// Returns iterator pointing to next element
			}
		}
	}

	/*	Discards the conic chords whose polynomial doesn't have any roots in [0, 1]
		Uses Sturm's algorithm for computing the number of roots of a polynomial
		Has crazy floating point problems
	*/
	void conicApproximationSturm(std::list<ConicChord>& chords, const ConicSection& conic,
								 Vector P10, Vector P20, Vector vP0,
						 		 Vector N10, Vector N20, Vector vN0,
								 Vector vL, Vector vV) const {
		std::list<ConicChord>::iterator it = chords.begin();
		
		while (it != chords.end()) {
			chordPolynomial(*it, m_eta, P10, P20, vP0, N10, N20, vN0, vL, vV);
			int sturmBound = sturmAlgorithm(it->polynomial, 0, 1);
			if (sturmBound == 0) {	// No roots, discard the chord
				it = chords.erase(it);
			} else if (sturmBound == 1) { // Single root, keep the chord
				++it;
			} else {	// More than one root over the chord [TODO]

				it = chords.erase(it);	// Returns iterator pointing to next element
			}
		}
	}

	/*	Poor man's algorithm for projecting a chord point onto the conic arc
		Uses the root t \in [0, 1] in a linear interpolation between t_0 and t_1
		(resp. the conic param of q0 and q1)
		It's a poor projection, see Eberly's modified method instead
	*/
	// [CHECK] Can the projected point not lie in the conic? (Ellipse case)
	Vector2 chordProjection(const ConicSection& conic, const ConicChord& chord,
							Float root) const {
		Float paramValue = (1.0 - root) * chord.q0.parameter.value +
								   root * chord.q1.parameter.value;
		ConicParam weightedParam = { .value = paramValue, .branch = chord.q0.parameter.branch };
		Vector3 projPos = conicParamToPos(conic, weightedParam);
		
		return Vector2(projPos[0], projPos[1]);
	}

	//-------------- End set of functions for conic sections ------------------

    //-------------- Begin set of functions for single scattering -------------
 	Spectrum attenuation(const Spectrum &muT, Float negDistance) const {
		Spectrum result(1.0);
		for (int c = 0; c < SPECTRUM_SAMPLES; ++c)
			if (m_sigmaT[c])
				result[c] = math::fastexp(muT[c] * negDistance);
		return result;
	}

	//-------------------------------------------------------------------------
	/* [MF2020] Modified aabbSegmentTest, but V is a single point
		instead of being in a segment. */
	bool aabbSpindleTest(const AABB &aabb, const Point &L,
						 const Point &V) const {
		// Is there a ray from V through any triangle
		// inside the AABB connecting to L?
		// Basic spindle test, but with cones.
		// Bounding sphere of the aabb:
		const BSphere triSphere = aabb.getBSphere();
		// Bounding cone for omega_L:
		Vector omegaL = L - triSphere.center;
		const Float domegaL = omegaL.length();
		omegaL /= domegaL;
		// Cone must be tangent to bounding sphere.
		if (domegaL < triSphere.radius)
			return true; // thetaL = M_PI. All tests will send true.

		const Float thetaL = math::safe_asin(triSphere.radius / domegaL);

		Vector omegaV = V - triSphere.center;
		const Float domegaV = omegaV.length();
		omegaV /= domegaV;

		if (domegaV < triSphere.radius)
			return true; // thetaH = M_PI. All tests will send true.

		const Float thetaV = math::safe_asin(triSphere.radius / domegaV);

		// Spindle test:
		// If omegaV is inside the cone of axis -omegaL, spindleAngle,
		// then there can be a spindle compatible intersection
		const Float spindleAngle =
			Float(0.5f * M_PI - math::safe_asin(m_invEta) + thetaL + thetaV);
		if (spindleAngle < M_PI) {
			const Float cosSpindle = std::cos(spindleAngle);
			Float cosVL = dot(omegaV, -omegaL);
			
			if (cosVL < cosSpindle)
				return false;
		}
		return true;
	}

	// [MF2020] [TODO] Code duplication in SpindleTest
	//------------------------------------------------------------------------
	/* [MF2020] Modified triangleSegmentTest, but V is a single point
		instead of being in a segment. */
	bool triangleSpindleTest(const Triangle &tri, const Point &L,
						 const Point &V, const Point *positions) const {
		// Is there a ray from V through triangle tri connecting to L?
		// Bounding sphere of the triangle:
		const BSphere triSphere = tri.getBSphere(positions);

		// Bounding cone for omega_L:
		Vector omegaL = L - triSphere.center;
		const Float domegaL = omegaL.length();
		omegaL /= domegaL;

		// Cone must be tangent to bounding sphere.
		if (domegaL < triSphere.radius)
			return true; // thetaL = M_PI. All tests will send true.

		const Float thetaL = math::safe_asin(triSphere.radius / domegaL);

		Vector omegaV = V - triSphere.center;
		const Float domegaV = omegaV.length();
		omegaV /= domegaV;

		if (domegaV < triSphere.radius)
			return true; // thetaH = M_PI. All tests will send true.

		const Float thetaV = math::safe_asin(triSphere.radius / domegaV);

		// Spindle test:
		// If omegaV is inside the cone of axis -omegaL, spindleAngle,
		// then there can be a spindle compatible intersection
		const Float spindleAngle =
			Float(0.5f * M_PI - math::safe_asin(m_invEta) + thetaL + thetaV);
		if (spindleAngle < M_PI) {
			const Float cosSpindle = std::cos(spindleAngle);
			Float cosVL = dot(omegaV, -omegaL);
			
			if (cosVL < cosSpindle)
				return false;
		}
		return true;
	}

	//------------------------------------------------------------------------
	/* [MF2020] Once the triangle passes the spindle test, we test for sidedness.
		If sidedness passes, then we compute the coplanarity conic and discard
		the triangle if the conic doesn't intersect it.
		
		If all of the tests pass, we compute the conic-triangle intersection points
		(there should be an even number of them unless the conic is tangent to the
		triangle at an intersection point). These intersections are paired such that
		a pair of points lie on the same chord of the conic section inside
		the triangle.

		[TODO] All functions that can fail should be checked (chordSplit among others)
	*/
	Spectrum testThisTriangle(const Triangle &tri, const Point &L,
							  const Point &V, const Vector &dInternal, Float dist,
							  const Point *positions, const Vector *normals,
							  const Spectrum &inputSpectrum, const Scene *scene,
							  Float time = 0.) const {
		const Point tP[3] = { positions[tri.idx[0]], positions[tri.idx[1]],
							  positions[tri.idx[2]] };
		const Vector P10 = tP[1] - tP[0];
		const Vector P20 = tP[2] - tP[0];

		Vector Ng = cross(P10, P20);
		const Float lNg = Ng.length();
		if (lNg < 1e-7f)
			return Spectrum(0.0f);
		Ng /= lNg;

		// Sidedness agreement, vs. geometric normal.
		// [TODO] Should also be checked for the shading normal
		const Vector PL0 = L - tP[0];
		if (dot(PL0, Ng) < 0)
			return Spectrum(0.0f);

		const Vector tN[3] = { normals[tri.idx[0]], normals[tri.idx[1]],
							   normals[tri.idx[2]] };
		const Vector N10 = tN[1] - tN[0];
		const Vector N20 = tN[2] - tN[0];

		// Point to Vector conversions for conic computation
		const Vector LV = L - V;
		const Vector vP0(tP[0]);
		const Vector vN0(tN[0]);
		const Vector vL(L);
		const Vector vV(V);

		// std::cout << "Triangle Parameters:" << std::endl;
		// std::cout << "P0:  " << tP[0].toString() << std::endl
		// 		  << "P1:  " << tP[1].toString() << std::endl
		// 		  << "P2:  " << tP[2].toString() << std::endl
		// 		  << "N0:  " << tN[0].toString() << std::endl
		// 		  << "N1:  " << tN[1].toString() << std::endl
		// 		  << "N2:  " << tN[2].toString() << std::endl
		// 		  << "L:   " << vL.toString() << std::endl
		// 		  << "V:   " << vV.toString() << std::endl
		// 		  << "eta: " << m_eta << std::endl;

		// Conic section coefficients in barycentric space
		const Float conicA = dot(N10, cross(P10, LV));
		const Float conicB = dot(N10, cross(P20, LV)) + dot(N20, cross(P10, LV));
		const Float conicC = dot(N20, cross(P20, LV));
		const Float conicD = dot(vN0, cross(P10, LV)) + dot(N10, cross(vP0, LV)) +
							 dot(N10, cross(vL, vV));
		const Float conicE = dot(vN0, cross(P20, LV)) + dot(N20, cross(vP0, LV)) +
							 dot(N20, cross(vL, vV));
		const Float conicF = dot(vN0, cross(vP0, LV)) + dot(vN0, cross(vL, vV));
		
		// std::cout << "Conic Parameters:" << std::endl;
		// std::cout << "A: " << conicA << std::endl
		// 		  << "B: " << conicB << std::endl
		// 		  << "C: " << conicC << std::endl
		// 		  << "D: " << conicD << std::endl
		// 		  << "E: " << conicE << std::endl
		// 		  << "F: " << conicF << std::endl;

		// 0, 1 --> b = 0 side (A)
		// 2, 3 --> a = 0 side (B)
		// 4, 5 --> a + b = 1 side AB)
		Float conicTriIsects[6];
		bool isUnitRoot[6];
		for (int i = 0; i < 6; ++i) {
			isUnitRoot[i] = false;
			/*	Initialize outside of [0, 1] so that even if there are no solutions
				for a given side of the triangle, the default values are not detected
				as roots in the [0, 1] interval. */
			conicTriIsects[i] = -1.0f;
		}

		// Find conic-triangle intersections in barycentric space
		solveQuadratic(conicA, conicD, conicF, conicTriIsects[0], conicTriIsects[1]);
		solveQuadratic(conicC, conicE, conicF, conicTriIsects[2], conicTriIsects[3]);
		solveQuadratic(conicA - conicB + conicC,
					   conicB - 2.0f * conicC + conicD - conicE,
					   conicC + conicE + conicF,
					   conicTriIsects[4], conicTriIsects[5]);
		
		// Find roots in [0, 1]
		// [TODO] Check that there is an even number of SINGLE roots (no multiplicity)
		bool hasUnitRoots = false;
		for (int i = 0; i < 6; ++i) {
			isUnitRoot[i] = (conicTriIsects[i] >= 0.0f) && (conicTriIsects[i] <= 1.0f);
			hasUnitRoots = hasUnitRoots || isUnitRoot[i];
		}

		// LAST TRIANGLE CULLING TEST
		// If the conic does not intersect the triangle, there cannot be any solutions
		// [TODO] Handle the case where the whole ellipse is in the triangle
		if (!hasUnitRoots) {
			return Spectrum(0.0f);
		}

		// START FINDING SOLUTIONS
		ConicSection conic ={ .coefficients = {conicA, conicB, conicC,
											   conicD, conicE, conicF}};

		// This function relies on numeric functions, it may fail
		if (!computeConicParameters(conic)) {
			Log(EError, "computeConicParameters() failed");
		}

		// Compute triangle-conic intersection points
		std::vector<IsectPoint> isectPoints;

		for (int i = 0; i < 6; ++i) {
			IsectPoint iPoint;
			ConicPoint cPoint;
			ETriangleSide triSide;
			// If the intersection is valid (root in [0,1])
			if (isUnitRoot[i]) {
				// Set side of triangle and compute barycentric coordinates of the point
				int side = i / 2;
				if (side == 0) {
					triSide = ETriangleSide::ESideA;
					cPoint.position = Vector(conicTriIsects[i], 0.0, 1.0);
				} else if (side == 1) {
					triSide = ETriangleSide::ESideB;
					cPoint.position = Vector(0.0, conicTriIsects[i], 1.0);
				} else if (side == 2) {
					triSide = ETriangleSide::ESideAB;
					cPoint.position = Vector(conicTriIsects[i], 1.0 - conicTriIsects[i], 1.0);
				}
				cPoint.parameter = conicPosToParam(conic, cPoint.position);
				iPoint.cPoint = cPoint;
				iPoint.side = triSide;
				isectPoints.push_back(iPoint);
			}
		}

		// Sort intersection points according to branch (Neg < Pos < Single) first
		// and parameter value second
		std::sort(isectPoints.begin(), isectPoints.end(), compIsectPts);
		int numIsects = isectPoints.size();
		ConicChord chord;
		std::list<ConicChord> conicApprox;

		// [CHECK] Proofread this section
		// [TODO] Deal with the case where the whole ellipse is contained in the triangle
		if (numIsects == 2) {
			// One chord
			// E/H: single pair of intersection points
			chord = {isectPoints[0].cPoint, isectPoints[1].cPoint};
			chordSplit(conic, chord);
			conicApprox.push_back(chord);
		} else if ((numIsects == 4) || (numIsects == 6)) {
			/*	Two or three chords
				[4]	E: pairs between adjacent intersection points not on the same triSide
					H: if single branch, same as for E; if two branches then pairs on same branches
				[6] E: same as 4 isects 
					H: single branch	*/
			if ((numIsects == 4) &&
				(isectPoints[0].cPoint.parameter.branch != isectPoints[3].cPoint.parameter.branch)) {
				// Two branches: [--++] configuration
				chord = {isectPoints[0].cPoint, isectPoints[1].cPoint};
				chordSplit(conic, chord);
				conicApprox.push_back(chord);
				chord = {isectPoints[2].cPoint, isectPoints[3].cPoint};
				chordSplit(conic, chord);
				conicApprox.push_back(chord);
			} else { // All other cases are single branch cases
				// Find the isectPoint sharing an edge with the next one in the list
				int sameSideIdx = 0;
				while (sameSideIdx < numIsects) {
					if (isectPoints[sameSideIdx].side == isectPoints[(sameSideIdx+1) % numIsects].side)
						break;
					++sameSideIdx;
				}

				int i = (sameSideIdx % 2 == 0) ? 1 : 0;
				while (i < numIsects) {
					chord = {isectPoints[i].cPoint, isectPoints[(i+1) % numIsects].cPoint};
					conicApprox.push_back(chord);
					i += 2;
				}	
			}
		} else {	// Should not happen (even with tangent points)
			Log(EError, "Odd number of intersection points");
		}

		/*	Subdivide the chords until the maximum distance between each chord and the
			conic falls below a user-specified threshold. Then discard every chord that
			doesn't contain a root, and subdivide those that contain multiple.	*/
		Float approxPrecision = 1e-10;
		conicApproximationSubdivide(conicApprox, conic, approxPrecision);
		conicApproximationDescartes(conicApprox, conic, P10, P20, vP0, N10, N20, vN0, vL, vV);

		Spectrum cfTri(0.0f);
		for (std::list<ConicChord>::iterator it = conicApprox.begin();
			 it != conicApprox.end(); ++it) {
			Float root = 0.0;
			if (!polynomialNewton(it->polynomial, root))	// Find the root on the chord
				Log(EError, "The root wasn't found in maxIter iterations");
			Vector2 conicRoot = chordProjection(conic, *it, root);	// Project the root on the conic

			Vector paramP(dist, conicRoot[0], conicRoot[1]);
			Float a11 = dot(P10, P10);
			Float a12 = dot(P10, P20);
			Float a22 = dot(P20, P20);
			const Float det = a11 * a22 - a12 * a12;
			if (det == 0) continue;	// [CHECK] comparison with 0??? Nicolas' code

			const Float invDet = 1.0f / det;
			a11 *= invDet;
			a12 *= invDet;
			a22 *= invDet;
	
			cfTri += contributionFromThatPoint(paramP, P10, P20, N10, N20, L, V, dInternal,
											   tP[0], tN, Ng, a11, a12, a22, inputSpectrum,
											   scene, time);
		}

		return cfTri;
	}

	//------------------------------------------------------------------------
	Spectrum contributionFromThatPoint(const Vector &paramP, const Vector &dPdu,
									   const Vector &dPdv, const Vector &dNsdu,
									   const Vector &dNsdv, const Point &L,
									   const Point &V, const Vector &dInternal,
									   const Point &P0, const Vector tN[3],
									   const Vector &Ng, Float a11, Float a12,
									   Float a22, const Spectrum &inputSpectrum,
									   const Scene *scene,
									   Float time = 0.0f) const {
		const Point Pc = P0 + paramP[1] * dPdu + paramP[2] * dPdv;
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

		// For debug/explanation only: result without the ray-differentials
		Float D;
		if (!m_distanceCorrection) {
			D = (domegaV + m_eta * domegaL) * (std::abs(cosThetaL/cosThetaV)*domegaV +
											 std::abs(cosThetaV/cosThetaL)*m_eta*domegaL);
			return result / D;
		}

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
		D = cross(L_p, L_s).length();

		//std::cout << "result / D:\t" << (result / D).toString() << std::endl;
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

		// [CHECK] the correctness of the computation of the contribution
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
		} else {
			// [MF2020] "Improved" algorithm start
			const Float sMax = 1 - exp(-(thickness / m_radius));
			const Float dSamples =
				m_singleScatterSamples
					? sMax / Float(m_singleScatterSamples)
					: sMax;

			// dRec.dist^2 is accounted for in inputSpectrum in testThisTriangle
			const Spectrum weight0 = (dSamples * m_radius) * m_sigmaS;

			const TriMesh *triMesh = static_cast<const TriMesh *>(its.shape);
			const Point *positions = triMesh->getVertexPositions();
			const Vector *normals = triMesh->getVertexNormals();

			size_t numTriangles = triMesh->getTriangleCount();
			bool *doneThisTriangleBefore = new bool[numTriangles];

			// For each sample along the ray
			for (int s = 0; s < m_singleScatterSamples; ++s) {
				Float sample = sampler->next1D() * sMax;
				sampler->advance();
				const Float dist = -math::fastlog(1 - sample) * m_radius;
				const Point V = its.p + dist * dInternal;
				if (EXPECT_NOT_TAKEN(dist > thickness))
					continue;
				
				for (size_t i = 0; i < numTriangles; i++) {
					doneThisTriangleBefore[i] = false;
				}

				/* Sampling weight, because sampling should have been
				 * chanel-dependant */
				Spectrum weight = weight0 * math::fastexp(m_invRadius * dist);

				// scan the kd-tree, cull all nodes not intersecting with segment,
				// keep going.
				struct spindleStackEntry {
					const ShapeKDTree::KDNode *node;
					AABB aabb;
				};

				spindleStackEntry stack[MTS_KD_MAXDEPTH];
				int stackPos = 0;
				const ShapeKDTree::KDNode *currNode = scene->getKDTree()->getRoot();
				AABB currNodeAABB = scene->getKDTree()->getTightAABB();

				while (currNode != NULL) {
					// if bounding sphere of AABB passes the spindle test
					if (aabbSpindleTest(currNodeAABB, L, V)) {
						// if node is a leaf
						if (currNode->isLeaf()) {
							const uint32_t primEnd = currNode->getPrimEnd();
							// for each primitive in the node
							for (uint32_t entry = currNode->getPrimStart();
								entry < primEnd; entry++) {
								uint32_t primIdx = scene->getKDTree()->getIndices()[entry];
								uint32_t shapeIdx = scene->getKDTree()->findShape(primIdx);
								if (its.shape ==
									scene->getKDTree()->getShapes()[shapeIdx]) {
									// it belongs to the original shape
									if (EXPECT_TAKEN(scene->getKDTree()->m_triangleFlag[shapeIdx])) {
										if (!doneThisTriangleBefore[primIdx]) {
											doneThisTriangleBefore[primIdx] = true;
											// if triangle passes spindle test
											if (triangleSpindleTest(triMesh->getTriangles()[primIdx],
													L, V, positions)) {
												result += weight *
													testThisTriangle(triMesh->getTriangles()[primIdx],
																	 L, V, dInternal, dist,
																	 positions, normals,
																	 value * (dRec.dist * dRec.dist),
																	 scene, its.time);
											}
										}
									}
								}
							}
							// Pop from the stack:
							if (stackPos > 0) {
								--stackPos;
								currNode = stack[stackPos].node;
								currNodeAABB = stack[stackPos].aabb;
							} else {
								break;
							}
						} else {	// node is not a leaf
							// add the two children to the stack, and iterate
							Float splitValue = currNode->getSplit();
							int axis = currNode->getAxis();
							stack[stackPos].node = currNode->getLeft();
							stack[stackPos].aabb = currNodeAABB;
							stack[stackPos].aabb.max[axis] = splitValue;
							currNode = currNode->getRight();
							currNodeAABB.min[axis] = splitValue;
							++stackPos;							
						}
					} else {	// spindle test fails
						// Pop from the stack:
						if (stackPos > 0) {
							--stackPos;
							currNode = stack[stackPos].node;
							currNodeAABB = stack[stackPos].aabb;
						} else {
							break;
						}
					}
				}
			}

			delete[] doneThisTriangleBefore;
		}
		return result;
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
	int m_singleScatterSamples;	// [MF2020]
};

MTS_IMPLEMENT_CLASS_S(SingleScatterOptimized, false, Subsurface)
MTS_EXPORT_PLUGIN(SingleScatterOptimized, "Single scattering model (new version)");
MTS_NAMESPACE_END
