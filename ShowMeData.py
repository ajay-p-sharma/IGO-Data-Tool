import os, sys, json, datetime
from flask import Flask, render_template, request, redirect, flash, url_for
from RevProxy import ReverseProxied
from werkzeug.utils import secure_filename
import requests
import os, json, re, yaml
from settings import APP_ROOT, APP_STATIC
app = Flask(__name__)
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.poolmanager import PoolManager
import ssl
import re
import time
import uwsgi, pickle
from operator import itemgetter

app.config['PROPAGATE_EXCEPTIONS'] = True
app.wsgi_app = ReverseProxied(app.wsgi_app)


class MyAdapter(HTTPAdapter):
    def init_poolmanager(self, connections, maxsize, block=False):
        self.poolmanager = PoolManager(num_pools=connections,
                maxsize=maxsize,
                block=block,
                ssl_version=ssl.PROTOCOL_SSLv23)
                #ssl_version=ssl.PROTOCOL_SSLv3)

s = requests.Session()
s.mount('https://', MyAdapter())


ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'])
UPLOAD_FOLDER = 'Uploads/'
recipeList = ["Agilent_MouseAllExonV1", "Agilent_v4_51MB_Human","AmpliconSeq","ChIPSeq","CHM","CRISPRSeq","CustomAmplificationPCR","CustomCapture",\
"HemePACT_v3", "HemePACT_v3+", "HemePACT_v4", "IDT_Exome_V1_IMPACT468","IDTCustomCaptur", "IMPACT410", "IMPACT410+","IMPACT468","IWG","IWGCustomCapture",\
"M-IMPACT_v1","NimblegenCustomCapture","R_Loop_DNA_Seq","RNASeq_PolyA","ShallowWGS","SMARTerAmpSeq","smRNASeq","WholeExomeSequencing","WholeGenomeBisulfateSequencing",\
"WholeGenomeSequencing"]
PATH_TO_STATS_FILE = "LimsStats.tsv"

@app.route('/')
def hello_world():
    datafile = PATH_TO_STATS_FILE
    fileNameSplit = datafile.split('/')
    filename = fileNameSplit[len(fileNameSplit)-1]

    if not os.path.isfile(datafile):
        flash('Data text file not found!')
        return redirect(url_for('fileNotFoundError'))

    if os.path.isfile(datafile):
        with open(datafile) as f:
            text = f.readlines()
            dataToShow = parseData(text)
            chartData = parse_chart_data(dataToShow)
            failure_data=get_library_failure_data_all(dataToShow)
            #impact_data = parseDataFilterAllByRecipe(text, "Impact410")
            #custom_capture_data = parseDataFilterAllByRecipe(text, "CustomCapture")
        return render_template('index.html', data = json.loads(dataToShow), data_to_chart = json.loads(chartData), failure_rate_data=json.loads(failure_data))


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/upload_file',  methods=['GET', 'POST'])
def upload_file():

    startDate = request.form['startDate']
    endDate = request.form['endDate']
    recipe = request.form['recipe']
    dataToShow = None
    datafile = PATH_TO_STATS_FILE
    fileNameSplit = datafile.split('/')
    filename = fileNameSplit[len(fileNameSplit) - 1]

    if not os.path.isfile(datafile):
        flash('Data text file not found!')
        return redirect(url_for('FileNotFound'))

    if os.path.isfile(datafile):

        with open(datafile) as f:
            text = f.readlines()

            if recipe != "" and startDate != "" and endDate != "":
                dataToShow = parseDataFilterByDateAndRecipe(text, startDate, endDate, recipe)
                chartData = parse_chart_data(dataToShow)
                failure_data=get_library_failure_data_all(dataToShow)
            if recipe == "All" and startDate != "" and endDate != "":
                dataToShow = parseDataFilterAllByDate(text, startDate, endDate)
                chartData = parse_chart_data(dataToShow)
                failure_data=get_library_failure_data_all(dataToShow)
            if recipe != "All" and recipe != "" and startDate == "" and endDate == "":
                dataToShow = parseDataFilterAllByRecipe(text, recipe)
                chartData = parse_chart_data(dataToShow)
                failure_data=get_library_failure_data_all(dataToShow)
            if recipe == "All" and startDate == "" and endDate == "":
                dataToShow = parseData(text)
                chartData = parse_chart_data(dataToShow)
                failure_data=get_library_failure_data_all(dataToShow)

            if getKeys(dataToShow) is True:
                return redirect(url_for('dataNotFound'))

            else:
                return render_template('LibraryPrepData.html', data=json.loads(dataToShow), data_to_chart=json.loads(chartData), failure_rate_data=json.loads(failure_data))



@app.route('/dataNotFound')
def dataNotFound():
    return render_template('DataNotFound.html', nodata="I cannot find the data for the selected Dates/Recipe.")


@app.route('/fileNotFound')
def fileNotFoundError():
    return render_template('FileNotFound.html', error="I cannot find the data file where I read from. Please make sure data file is present.")


def parseData(data):
    dictMonths = {'January': "01", 'February': "02", 'March': "03", 'April': "04", 'May': "05", 'June': "06", 'July': "07", 'August': "08",
                  'September': "09", 'October': "10", 'November': "11", 'December': "12"}
    month = ""
    year = 0.0
    initialInput = 0.0
    inputUnits = ""
    recipe = ""
    preservation = ""
    volume = ""
    concentration = ""
    qcStatus = ""
    libraryYield = 0.0
    finalData = {}
    count = 0
    datad = {}


    for line in data:
        values = line.split('\t')

        if values[2].find("DNALibraryPrepProtocol1") !=-1 and values[3]=="OtherSampleId":

            initialInput = 0.0
            inputUnits = ""
            recipe = ""
            volume = 0.0
            qcStatus = ""
            libraryYield = 0.0


            if datad:

                monthYearKey = str(year)+ '/' + str(dictMonths[month])

                finalData.setdefault(monthYearKey, {})
                finalData[monthYearKey].setdefault('libraryData', []).append(datad)
                count = count + 1

                if 'minInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minInput'] = datad['input']

                if float(datad['input']) < float(finalData[monthYearKey]['minInput']):
                    finalData[monthYearKey]['minInput'] = datad['input']

                if 'maxInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if float(datad['input']) > float(finalData[monthYearKey]['maxInput']):
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if 'sumInput' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumInput'] = float(datad['input'])

                else:

                    finalData[monthYearKey]['sumInput'] += datad['input']

                if 'averageInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageInput'] = datad['input']

                else:

                    finalData[monthYearKey]['averageInput'] = finalData[monthYearKey]['sumInput']/len(finalData[monthYearKey]['libraryData'])

                if 'minLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if float(datad['libraryYield']) < float(finalData[monthYearKey]['minLibYield']):
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if 'maxLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxLibYield'] = libraryYield

                if float(datad['libraryYield']) > float(finalData[monthYearKey]['maxLibYield']):
                    finalData[monthYearKey]['maxLibYield'] = datad['libraryYield']

                if 'sumLibYield' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumLibYield'] = float(datad['libraryYield'])

                else:

                    finalData[monthYearKey]['sumLibYield'] += float(datad['libraryYield'])

                if 'averageLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageLibYield'] = libraryYield

                else:

                    finalData[monthYearKey]['averageLibYield'] = float(finalData[monthYearKey]['sumLibYield'])/len(finalData[monthYearKey]['libraryData'])

                if datad['recipe'] != "" and datad['recipe'] not in recipeList:
                    recipeList.append(datad['recipe'])

                datad = {}

            date = datetime.datetime.fromtimestamp(float(values[0])/1000).strftime('%Y-%m-%d %H:%M:%S')
            month=datetime.datetime.fromtimestamp(float(values[0])/1000).strftime('%B')
            year=datetime.datetime.fromtimestamp(float(values[0])/1000).strftime('%Y')

        elif values[2].find("QCProtocol") !=-1 and values[3] == "Concentration":
            sample = values[2].split(':')[1]

        if values[2].find("DNALibraryPrepProtocol1") !=-1 and values[3] == "TargetMassAliq1":

            initialInput = values[4].strip()

            if initialInput == 'null':
                initialInput = 0.01

            inputUnits = "ng"

        if values[2].find("Sample") != -1 and values[3] == "Recipe":
            recipe = values[4]

        if values[2].find("Sample") != -1 and values[3] == "Preservation":
            preservation = values[4]

        if values[2].find("QCProtocol") !=-1 and values[3] == "Concentration":
            concentration = values[4]
            sample = values[2].split(':')[1]

        if values[2].find("QCProtocol") !=-1 and values[3] == "ConcentrationUnits":
            concentrationUnits = values[4]

        if values[2].find("QCProtocol") != -1 and values[3] == "Volume":
            volume = values[4]
            if concentration.strip()=='null' or volume.strip() =='null':
                concentration=0.0
                volume = 0.0
               
                libraryYield = float(concentration)*float(volume)
            else:
                libraryYield = float(concentration) * float(volume)

        if values[2].find("QCProtocol") !=-1 and values[3] == "FinalQcStatus":
            qcStatus=values[4]

        if initialInput >= 0.0 and recipe != "" and qcStatus != "" and volume >= 0.0 and libraryYield >= 0.0:

            datad['sample'] = sample
            datad['month'] = month.strip()
            datad['year'] = year.strip()
            datad['input'] = float(initialInput)
            datad['recipe'] = recipe.strip()
            datad['preservation'] = preservation.strip()
            datad['concentration'] = concentration
            datad['concentrationUnits'] = concentrationUnits.strip()
            datad['volume'] = volume
            datad['libraryYield'] = float(libraryYield)
            datad['qcStatus'] = qcStatus.strip()
            datad['inputUnits'] = inputUnits.strip()
            datad['count'] = int(count)

    finalData['recipeList'] = recipeList
   
    return json.dumps(finalData, indent=4, sort_keys=True)


def parseDataFilterByDateAndRecipe(data, startDate, endDate, Recipe):
    dictMonths = {'January': "01", 'February': "02", 'March': "03", 'April': "04", 'May': "05", 'June': "06",
                  'July': "07", 'August': "08", 'September': "09", 'October': "10", 'November': "11", 'December': "12"}

    startMonth = int(startDate.split('/')[0])
    startYear = int(startDate.split('/')[1])
    endMonth = int(endDate.split('/')[0])
    endYear = int(endDate.split('/')[1])
    recipeVal = Recipe
    month = ""
    year = 0.0
    initialInput = 0.0
    inputUnits = ""
    recipe = ""
    preservation = ""
    volume = 0.0
    qcStatus = ""
    libraryYield = 0.0
    finalData = {}
    count = 0
    datad = {}

    for line in data:
        values = line.split('\t')

        if values[2].find("DNALibraryPrepProtocol1") != -1 and values[3] == "OtherSampleId":

            initialInput = 0.0
            inputUnits = ""
            recipe = ""
            volume = 0.0
            qcStatus = ""
            libraryYield = 0.0


            if datad:
                monthYearKey=str(year)+ '/' + str(dictMonths[month])
                finalData.setdefault(monthYearKey, {})
                finalData[monthYearKey].setdefault('libraryData', []).append(datad)
                count = count + 1

                if 'minInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minInput'] = datad['input']

                if float(datad['input']) < float(finalData[monthYearKey]['minInput']):
                    finalData[monthYearKey]['minInput'] = datad['input']

                if 'maxInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if float(datad['input']) > float(finalData[monthYearKey]['maxInput']):
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if 'sumInput' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumInput'] = float(datad['input'])

                else:

                    finalData[monthYearKey]['sumInput'] += datad['input']

                if 'averageInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageInput'] = datad['input']

                else:

                    finalData[monthYearKey]['averageInput'] = finalData[monthYearKey]['sumInput'] / len(
                        finalData[monthYearKey]['libraryData'])

                if 'minLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if float(datad['libraryYield']) < float(finalData[monthYearKey]['minLibYield']):
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if 'maxLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxLibYield'] = libraryYield

                if float(datad['libraryYield']) > float(finalData[monthYearKey]['maxLibYield']):
                    finalData[monthYearKey]['maxLibYield'] = datad['libraryYield']

                if 'sumLibYield' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumLibYield'] = float(datad['libraryYield'])

                else:

                    finalData[monthYearKey]['sumLibYield'] += float(datad['libraryYield'])

                if 'averageLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageLibYield'] = libraryYield

                else:
                    finalData[monthYearKey]['averageLibYield'] = float(finalData[monthYearKey]['sumLibYield']) / len(
                    finalData[monthYearKey]['libraryData'])

                datad = {}
            date = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%Y-%m-%d %H:%M:%S')
            month = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%B')
            year = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%Y')
            year=int(year)

        if values[2].find("DNALibraryPrepProtocol1") != -1 and values[3] == "TargetMassAliq1":
            initialInput = values[4].strip()

            if initialInput == 'null':
                initialInput = 0.01
            inputUnits = "ng"

        if values[2].find("Sample") != -1 and values[3] == "Recipe":
            recipe = values[4]

        if values[2].find("Sample") != -1 and values[3] == "Preservation":
            preservation = values[4]

        if values[2].find("QCProtocol") != -1 and values[3] == "Concentration" and values[4] != "null":
            concentration = values[4]
            sample = values[2].split(':')[1]

        if values[2].find("QCProtocol") != -1 and values[3] == "ConcentrationUnits":
            concentrationUnits = values[4]

        if values[2].find("QCProtocol") != -1 and values[3] == "Volume":
            volume = values[4]
            if concentration.strip()=='null' or volume.strip() =='null':
                concentration=0.0
                volume = 0.0
               
                libraryYield = float(concentration)*float(volume)
            else:
                libraryYield = float(concentration) * float(volume)

        if values[2].find("QCProtocol") != -1 and values[3] == "FinalQcStatus":
            qcStatus = values[4]

        if initialInput >= 0.0 and recipe != "" and str(recipe.lower().strip()) == str(recipeVal.lower().strip()) and qcStatus != "" and volume >= 0.0 \
                    and libraryYield >= 0.0 and int(dictMonths[month]) >= startMonth and int(dictMonths[month]) <= endMonth and year>=startYear and year<=endYear:
            datad['sample'] = sample
            datad['month'] = month.strip()
            datad['year'] = year
            datad['input'] = float(initialInput)
            datad['recipe'] = recipe.strip()
            datad['preservation'] = preservation.strip()
            datad['concentration'] = float(concentration)
            datad['concentrationUnits'] = concentrationUnits.strip()
            datad['volume'] = float(volume)
            datad['libraryYield'] = float(libraryYield)
            datad['qcStatus'] = qcStatus.strip()
            datad['inputUnits'] = inputUnits.strip()
            datad['count'] = int(count)


    finalData['recipeList'] = recipeList
    return json.dumps(finalData, indent=4, sort_keys=True)

def parseDataFilterAllByDate(data, startDate, endDate):
    dictMonths = {'January': "01", 'February': "02", 'March': "03", 'April': "04", 'May': "05", 'June': "06",
                  'July': "07", 'August': "08", 'September': "09", 'October': "10", 'November': "11", 'December': "12"}
    startMonth= int(startDate.split('/')[0])
    startYear=int(startDate.split('/')[1])
    endMonth = int(endDate.split('/')[0])
    endYear = int(endDate.split('/')[1])
    month = ""
    year = 0.0
    initialInput = 0.0
    inputUnits = ""
    recipe = ""
    preservation = ""
    volume = 0.0
    qcStatus = ""
    libraryYield = 0.0
    finalData = {}
    count = 0
    datad = {}

    for line in data:
        values = line.split('\t')

        if (values[2].find("DNALibraryPrepProtocol1") != -1 and values[3] == "OtherSampleId"):

            initialInput = 0.0
            inputUnits = ""
            recipe = ""
            volume = 0.0
            qcStatus = ""
            libraryYield = 0.0


            if datad:
                monthYearKey=str(year)+ '/' + str(dictMonths[month])
                finalData.setdefault(monthYearKey, {})
                finalData[monthYearKey].setdefault('libraryData', []).append(datad)
                count = count + 1

                if 'minInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minInput'] = datad['input']

                if float(datad['input']) < float(finalData[monthYearKey]['minInput']):
                    finalData[monthYearKey]['minInput'] = datad['input']

                if 'maxInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if float(datad['input']) > float(finalData[monthYearKey]['maxInput']):
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if 'sumInput' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumInput'] = float(datad['input'])

                else:

                    finalData[monthYearKey]['sumInput'] += datad['input']

                if 'averageInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageInput'] = datad['input']

                else:

                    finalData[monthYearKey]['averageInput'] = finalData[monthYearKey]['sumInput'] / len(
                        finalData[monthYearKey]['libraryData'])

                if 'minLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if float(datad['libraryYield']) < float(finalData[monthYearKey]['minLibYield']):
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if 'maxLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxLibYield'] = libraryYield

                if float(datad['libraryYield']) > float(finalData[monthYearKey]['maxLibYield']):
                    finalData[monthYearKey]['maxLibYield'] = datad['libraryYield']

                if 'sumLibYield' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumLibYield'] = float(datad['libraryYield'])

                else:

                    finalData[monthYearKey]['sumLibYield'] += float(datad['libraryYield'])

                if 'averageLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageLibYield'] = libraryYield

                else:
                    finalData[monthYearKey]['averageLibYield'] = float(finalData[monthYearKey]['sumLibYield']) / len(
                    finalData[monthYearKey]['libraryData'])

                datad = {}
            date = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%Y-%m-%d %H:%M:%S')
            month = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%B')
            year = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%Y')
            year=int(year)

        if values[2].find("DNALibraryPrepProtocol1") != -1 and values[3] == "TargetMassAliq1":
            initialInput = values[4].strip()
            if initialInput == 'null':
                initialInput = 0.01
            inputUnits = "ng"

        if values[2].find("Sample") != -1 and values[3] == "Recipe":
            recipe = values[4]

        if values[2].find("Sample") != -1 and values[3] == "Preservation":
            preservation = values[4]

        if values[2].find("QCProtocol") != -1 and values[3] == "Concentration":
            concentration = values[4]
            sample = values[2].split(':')[1]

        if values[2].find("QCProtocol") != -1 and values[3] == "ConcentrationUnits":
            concentrationUnits = values[4]

        if values[2].find("QCProtocol") != -1 and values[3] == "Volume":
            volume = values[4]
            if concentration.strip()=='null' or volume.strip() =='null':
                concentration=0.0
                volume = 0.0
                
                libraryYield = float(concentration)*float(volume)
            else:
                libraryYield = float(concentration) * float(volume)

        if values[2].find("QCProtocol") != -1 and values[3] == "FinalQcStatus":
            qcStatus = values[4]

        if initialInput >= 0.0 and recipe != "" and qcStatus != "" and volume >= 0.0 and libraryYield >= 0.0 and int(dictMonths[month]) >= startMonth and int(dictMonths[month]) <= endMonth and year >= startYear and year <= endYear:
            datad['sample'] = sample
            datad['month'] = month.strip()
            datad['year'] = year
            datad['input'] = float(initialInput)
            datad['recipe'] = recipe.strip()
            datad['preservation'] = preservation.strip()
            datad['concentration'] = float(concentration)
            datad['concentrationUnits'] = concentrationUnits.strip()
            datad['volume'] = float(volume)
            datad['libraryYield'] = float(libraryYield)
            datad['qcStatus'] = qcStatus.strip()
            datad['inputUnits'] = inputUnits.strip()
            datad['count'] = int(count)


    finalData['recipeList'] = recipeList

    

    return json.dumps(finalData, indent=4, sort_keys=True)


def parseDataFilterAllByRecipe(data, Recipe):
    dictMonths = {'January': "01", 'February': "02", 'March': "03", 'April': "04", 'May': "05", 'June': "06",
                  'July': "07", 'August': "08", 'September': "09", 'October': "10", 'November': "11", 'December': "12"}
    recipeVal = Recipe
    month = ""
    year = 0.0
    initialInput = 0.0
    inputUnits = ""
    recipe = ""
    preservation = ""
    volume = 0.0
    qcStatus = ""
    libraryYield = 0.0
    finalData = {}
    count = 0
    datad = {}

    for line in data:
        values = line.split('\t')

        if values[2].find("DNALibraryPrepProtocol1") != -1 and values[3] == "OtherSampleId":

            initialInput = 0.0
            inputUnits = ""
            recipe = ""
            volume = 0.0
            qcStatus = ""
            libraryYield = 0.0


            if datad:
                monthYearKey=str(year)+ '/' + str(dictMonths[month])
                finalData.setdefault(monthYearKey, {})
                finalData[monthYearKey].setdefault('libraryData', []).append(datad)
                count = count + 1

                if 'minInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minInput'] = datad['input']

                if float(datad['input']) < float(finalData[monthYearKey]['minInput']):
                    finalData[monthYearKey]['minInput'] = datad['input']

                if 'maxInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if float(datad['input']) > float(finalData[monthYearKey]['maxInput']):
                    finalData[monthYearKey]['maxInput'] = datad['input']

                if 'sumInput' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumInput'] = float(datad['input'])

                else:

                    finalData[monthYearKey]['sumInput'] += datad['input']

                if 'averageInput' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageInput'] = datad['input']

                else:

                    finalData[monthYearKey]['averageInput'] = finalData[monthYearKey]['sumInput'] / len(
                        finalData[monthYearKey]['libraryData'])

                if 'minLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if float(datad['libraryYield']) < float(finalData[monthYearKey]['minLibYield']):
                    finalData[monthYearKey]['minLibYield'] = datad['libraryYield']

                if 'maxLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['maxLibYield'] = libraryYield

                if float(datad['libraryYield']) > float(finalData[monthYearKey]['maxLibYield']):
                    finalData[monthYearKey]['maxLibYield'] = datad['libraryYield']

                if 'sumLibYield' not in finalData[monthYearKey]:

                    finalData[monthYearKey]['sumLibYield'] = float(datad['libraryYield'])

                else:

                    finalData[monthYearKey]['sumLibYield'] += float(datad['libraryYield'])

                if 'averageLibYield' not in finalData[monthYearKey]:
                    finalData[monthYearKey]['averageLibYield'] = libraryYield

                else:
                    finalData[monthYearKey]['averageLibYield'] = float(finalData[monthYearKey]['sumLibYield']) / len(
                    finalData[monthYearKey]['libraryData'])

                datad = {}
            date = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%Y-%m-%d %H:%M:%S')
            month = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%B')
            year = datetime.datetime.fromtimestamp(float(values[0]) / 1000).strftime('%Y')

        if values[2].find("DNALibraryPrepProtocol1") != -1 and values[3] == "TargetMassAliq1":
            initialInput = values[4].strip()
            if initialInput == 'null':
                initialInput = 0.01
            inputUnits = "ng"

        if values[2].find("Sample") != -1 and values[3] == "Recipe":
            recipe = values[4]

        if values[2].find("Sample") != -1 and values[3] == "Preservation":
            preservation = values[4]

        if values[2].find("QCProtocol") != -1 and values[3] == "Concentration":
            concentration = values[4]
            sample = values[2].split(':')[1]

        if values[2].find("QCProtocol") != -1 and values[3] == "ConcentrationUnits":
            concentrationUnits = values[4]

        if values[2].find("QCProtocol") != -1 and values[3] == "Volume":
            volume = values[4]
            if concentration.strip()=='null' or volume.strip() =='null':
                concentration=0.0
                volume = 0.0
               
                libraryYield = float(concentration)*float(volume)
            else:
                libraryYield = float(concentration) * float(volume)

        if values[2].find("QCProtocol") != -1 and values[3] == "FinalQcStatus":
            qcStatus = values[4]

        if initialInput >= 0.0 and recipe != "" and str(recipe.lower().strip()) == str(recipeVal.lower().strip()) and qcStatus != "" and volume >= 0.0 and libraryYield >= 0.0:
            datad['sample'] = sample
            datad['month'] = month.strip()
            datad['year'] = year.strip()
            datad['input'] = float(initialInput)
            datad['recipe'] = recipe.strip()
            datad['preservation'] = preservation.strip()
            datad['concentration'] = float(concentration)
            datad['concentrationUnits'] = concentrationUnits.strip()
            datad['volume'] = float(volume)
            datad['libraryYield'] = float(libraryYield)
            datad['qcStatus'] = qcStatus.strip()
            datad['inputUnits'] = inputUnits.strip()
            datad['count'] = int(count)

    finalData['recipeList'] = recipeList

    
    return json.dumps(finalData, indent=4, sort_keys=True)


def parse_chart_data(data):
    labels=["month/year", "Average Input", "Average Yield"]
    temp_data=[]
    chartData=[]
    chartData.append(labels)
    data= json.loads(data)

    for key, values in sorted(data.iteritems()):

        if key != "recipeList":
            temp_data.append(key)
            temp_data.append(data[key]['averageInput'])
            temp_data.append(data[key]['averageLibYield'])
            chartData.append(temp_data)
            temp_data=[]
   

    return json.dumps(chartData, indent=4, sort_keys=True)


def getKeys(data):
    toTestData = json.loads(data)
   
    keysArray=[]
    for key, values in toTestData.items():
        keysArray.append(key)

    return len(keysArray) <= 1

def get_library_failure_data_all(data):
    data=json.loads(data)
    failed_count=0
    total_libraries_prepared=0
    failure_aggregate_data={}
    for key,value in data.items():
        if key =='recipeList':
            pass
        if 'libraryData' in value:
            for items in value['libraryData']:
                if items['qcStatus']== "Failed" or items['qcStatus']=="Failed - Reprocess":
                    failed_count = failed_count + 1

            failure_aggregate_data[key]=(float(failed_count)/float(len(value['libraryData'])))*100
            failed_count=0 

   
    return json.dumps(failure_aggregate_data, indent=4, sort_keys=True)
    
    
        
if __name__ == '__main__':
    app.secret_key = 'super secret key'
    
    app.run(debug = False)
