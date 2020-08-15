import subprocess

inputPath = 'scenes/input/bumpysphere/'
outputPath = 'scenes/output/bumpysphere/'
logPath = 'scenes/logs/bumpysphere/'

ssPath = 'singlescatter_path/'
ssOptPath = 'singlescatteropt_path/'
inputFile = 'bumpysphere.xml'

# In order: fastSingleScatter, distanceCorrection, maxDepth, sampleCount
testValues = [['false', 'true', '5', '1']]

def generateFlags(paramList):
    flags = []
    flags.append('-v')
    
    flags.append('-DfastSS=' + paramList[0])
    flags.append('-DdistC=' + paramList[1])    
    flags.append('-Ddepth=' + paramList[2])
    flags.append('-Dsamples=' + paramList[3])

    return flags

def shortBool(boolStr):
    if boolStr == 'true':
        return 'T'
    return 'F'

def generateOutputName(paramList):
    return 'bumpysphere' + '_fss' + shortBool(paramList[0]) + '_dcr' + shortBool(paramList[1]) + '_dpt' + paramList[2] + '_smp' + paramList[3]

def runMitsuba(paramListList):
    #subprocess.run(['source', 'setpath.sh'])
    
    inputSS = inputPath + ssPath + inputFile
    inputSSO = inputPath + ssOptPath + inputFile
    
    for pList in paramListList:
        commandSS = []
        commandSSO = []
        
        commandSS.append('mitsuba')
        commandSSO.append('mitsuba')

        flags = generateFlags(pList)
        commandSS += flags
        commandSSO += flags

        commandSS.append(inputSS)
        commandSSO.append(inputSSO)

        commandSS.append('-o')
        commandSSO.append('-o')

        outputFile = generateOutputName(pList)
        commandSS.append(outputPath + ssPath + outputFile)
        commandSSO.append(outputPath + ssOptPath + outputFile)

        print(" ".join(commandSS))
#        resultSS = subprocess.run(commandSS, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        print('\n')
        
        print(" ".join(commandSSO))
#        resultSSO = subprocess.run(commandSSO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)

#        with open(logPath + ssPath + outputFile + '.log', 'w') as ssFile:
#            ssFile.write(resultSS.stdout)
#        with open(logPath + ssOptPath + outputFile + '.log', 'w') as ssOptFile:
#            ssOptFile.write(resultSSO.stdout)

runMitsuba(testValues)
