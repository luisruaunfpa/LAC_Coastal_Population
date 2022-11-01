"""
Model exported as python.
Name : Lines to buffer process
Group : LAC Coastal population
With QGIS : 32210
"""

from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingMultiStepFeedback
from qgis.core import QgsProcessingParameterVectorLayer
from qgis.core import QgsProcessingParameterNumber
from qgis.core import QgsProcessingParameterBoolean
from qgis.core import QgsProcessingParameterFeatureSink
import processing


class LinesToBufferProcess(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterVectorLayer('adminbound', 'cline_raw', types=[QgsProcessing.TypeVectorLine], defaultValue=None))
        self.addParameter(QgsProcessingParameterNumber('splitlineparameter', 'Split line parameter', type=QgsProcessingParameterNumber.Double, defaultValue=100))
        self.addParameter(QgsProcessingParameterNumber('10kmbuffdist', '10km_buff_dist', type=QgsProcessingParameterNumber.Double, defaultValue=10000))
        self.addParameter(QgsProcessingParameterNumber('1kmbuffdist', '1km_buff_dist', type=QgsProcessingParameterNumber.Double, defaultValue=1000))
        self.addParameter(QgsProcessingParameterNumber('5kmbuffdist', '5km_buff_dist', type=QgsProcessingParameterNumber.Double, defaultValue=5000))
        self.addParameter(QgsProcessingParameterBoolean('buffdissolvecondition', 'buff dissolve condition', defaultValue=True))
        self.addParameter(QgsProcessingParameterNumber('segmentsbuff', 'segments buff ', type=QgsProcessingParameterNumber.Double, defaultValue=20))
        self.addParameter(QgsProcessingParameterFeatureSink('Buff1_54034', 'buff1_54034', type=QgsProcessing.TypeVectorPolygon, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Buff5_54034', 'buff5_54034', type=QgsProcessing.TypeVectorPolygon, createByDefault=True, supportsAppend=True, defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('Buff10_54034', 'buff10_54034', type=QgsProcessing.TypeVectorPolygon, createByDefault=True, supportsAppend=True, defaultValue=None))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(5, model_feedback)
        results = {}
        outputs = {}

        # Multipart to singleparts
        alg_params = {
            'INPUT': parameters['adminbound'],
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['MultipartToSingleparts'] = processing.run('native:multiparttosingleparts', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}

        # Split lines by maximum length
        # Split lines every 1km to avoid gap in buffers artifacts
        alg_params = {
            'INPUT': outputs['MultipartToSingleparts']['OUTPUT'],
            'LENGTH': parameters['splitlineparameter'],
            'OUTPUT': QgsProcessing.TEMPORARY_OUTPUT
        }
        outputs['SplitLinesByMaximumLength'] = processing.run('native:splitlinesbylength', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(2)
        if feedback.isCanceled():
            return {}

        # Buffer 5km
        alg_params = {
            'DISSOLVE': parameters['buffdissolvecondition'],
            'DISTANCE': parameters['10kmbuffdist'],
            'END_CAP_STYLE': 0,  # Round
            'INPUT': outputs['SplitLinesByMaximumLength']['OUTPUT'],
            'JOIN_STYLE': 0,  # Round
            'MITER_LIMIT': 2,
            'SEGMENTS': parameters['segmentsbuff'],
            'OUTPUT': parameters['Buff10_54034']
        }
        outputs['Buffer5km'] = processing.run('native:buffer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Buff10_54034'] = outputs['Buffer5km']['OUTPUT']

        feedback.setCurrentStep(3)
        if feedback.isCanceled():
            return {}

        # Buffer 1km
        alg_params = {
            'DISSOLVE': parameters['buffdissolvecondition'],
            'DISTANCE': parameters['1kmbuffdist'],
            'END_CAP_STYLE': 0,  # Round
            'INPUT': outputs['SplitLinesByMaximumLength']['OUTPUT'],
            'JOIN_STYLE': 0,  # Round
            'MITER_LIMIT': 2,
            'SEGMENTS': parameters['segmentsbuff'],
            'OUTPUT': parameters['Buff1_54034']
        }
        outputs['Buffer1km'] = processing.run('native:buffer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Buff1_54034'] = outputs['Buffer1km']['OUTPUT']

        feedback.setCurrentStep(4)
        if feedback.isCanceled():
            return {}

        # Buffer 5km
        alg_params = {
            'DISSOLVE': parameters['buffdissolvecondition'],
            'DISTANCE': parameters['5kmbuffdist'],
            'END_CAP_STYLE': 0,  # Round
            'INPUT': outputs['SplitLinesByMaximumLength']['OUTPUT'],
            'JOIN_STYLE': 0,  # Round
            'MITER_LIMIT': 2,
            'SEGMENTS': parameters['segmentsbuff'],
            'OUTPUT': parameters['Buff5_54034']
        }
        outputs['Buffer5km'] = processing.run('native:buffer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Buff5_54034'] = outputs['Buffer5km']['OUTPUT']
        return results

    def name(self):
        return 'Lines to buffer process'

    def displayName(self):
        return 'Lines to buffer process'

    def group(self):
        return 'LAC Coastal population'

    def groupId(self):
        return 'LAC Coastal population'

    def createInstance(self):
        return LinesToBufferProcess()
