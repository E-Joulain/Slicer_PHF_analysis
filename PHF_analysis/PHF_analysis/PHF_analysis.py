import logging
import os
from typing import Annotated, Optional
from scipy.optimize import minimize

import vtk
import qt
import numpy as np
import math
import ctk
import csv
import ast

import slicer
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from slicer.parameterNodeWrapper import (
    parameterNodeWrapper,
    WithinRange,)
from slicer import vtkMRMLScalarVolumeNode

try:
    import pandas as pd
except ImportError:
    slicer.util.infoDisplay("Installing pandas. Please wait...")
    with slicer.util.displayPythonShell():
        slicer.util.pip_install('pandas')
    import pandas as pd


#
# PHF_analysis
#


class PHF_analysis(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("PHF_analysis")  # TODO: make this more human readable by adding spaces
        # TODO: set categories (folders where the module shows up in the module selector)
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "Costume_model")]
        self.parent.dependencies = []  # TODO: add here list of module names that this module requires
        self.parent.contributors = ["Elise Joulain (UNIGE)"]  # TODO: replace with "Firstname Lastname (Organization)"
        # TODO: update with short description of the module and a link to online module documentation
        # _() function marks text as translatable to other languages
        self.parent.helpText = _("""
This is an example of scripted loadable module bundled in an extension.
See more information in <a href="https://github.com/organization/projectname#PHF_analysis">module documentation</a>.
""")
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = _("""
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""")

        # Additional initialization step after application startup is complete
        slicer.app.connect("startupCompleted()", registerSampleData)


#
# Register sample data sets in Sample Data module
#


def registerSampleData():
    """Add data sets to Sample Data module."""
    # It is always recommended to provide sample data for users to make it easy to try the module,
    # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

    import SampleData

    iconsPath = os.path.join(os.path.dirname(__file__), "Resources/Icons")

    # To ensure that the source code repository remains small (can be downloaded and installed quickly)
    # it is recommended to store data sets that are larger than a few MB in a Github release.

    # PHF_analysis1
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category="PHF_analysis",
        sampleName="PHF_analysis1",
        # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
        # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
        thumbnailFileName=os.path.join(iconsPath, "PHF_analysis1.png"),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
        fileNames="PHF_analysis1.nrrd",
        # Checksum to ensure file integrity. Can be computed by this command:
        #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
        checksums="SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
        # This node name will be used when the data set is loaded
        nodeNames="PHF_analysis1",
    )

    # PHF_analysis2
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category="PHF_analysis",
        sampleName="PHF_analysis2",
        thumbnailFileName=os.path.join(iconsPath, "PHF_analysis2.png"),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
        fileNames="PHF_analysis2.nrrd",
        checksums="SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
        # This node name will be used when the data set is loaded
        nodeNames="PHF_analysis2",
    )







class PHF_analysisWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self.inputVolumeSelector = None
        self.markupNode = None

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        self.fiducialLabels = {
            "P1": "Distal point, lateral internal cortex",
            "P2": "Proximal point, lateral internal cortex",
            "P3": "Proximal point, medial internal cortex",
            "P4": "Distal point, medial internal cortex",
            "P5": "Medial diaphyseal reduction criterion",
            "P6": "Medial head fragment reduction criterion",
            "P7": "Medial bone-cartilage junction",
            "P8": "Lateral bone-cartilage junction",
            "P9": "Proximal greater tuberosity reduction criterion",
            "P10": "Greater tuberosity apex",
            "P11": "Head sphericity point 1",
            "P12": "Head sphericity point 2",
            "P13": "Head sphericity point 3",
        }

        formLayout = qt.QFormLayout()
        self.layout.addLayout(formLayout)
        # --- Input Volume Selector ---
        self.inputVolumeSelector = slicer.qMRMLNodeComboBox()
        self.inputVolumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.inputVolumeSelector.setMRMLScene(slicer.mrmlScene)
        self.inputVolumeSelector.currentNodeChanged.connect(self.onInputVolumeChanged)
        self.layout.addWidget(qt.QLabel("Choose the volume on which you want to work on:"))
        self.layout.addWidget(self.inputVolumeSelector)
        self.inputVolumeSelector.setStyleSheet("""
                    QComboBox {
                        background-color: #F44336;
                        color: #00363A;
                        font-weight: bold;
                        border: 1px solid #00796B;
                        padding: 3px;
                        border-radius: 4px;
                    }
                """)

        # --- Side Selection ---
        self.selectionCollapsibleButton = ctk.ctkCollapsibleButton()
        self.selectionCollapsibleButton.text = "Enter the side of the inspected shoulder"
        self.layout.addWidget(self.selectionCollapsibleButton)
        self.selectionFormLayout = qt.QFormLayout(self.selectionCollapsibleButton)
        self.sideSelector = qt.QComboBox()
        self.sideSelector.addItem("Right")
        self.sideSelector.addItem("Left")
        self.selectionFormLayout.addRow("Select Side:", self.sideSelector)

        # --- Patient ID Input ---
        self.patientIdLineEdit = qt.QLineEdit()
        self.patientIdLineEdit.setPlaceholderText("Enter patient ID")
        formLayout.addRow("Patient ID:", self.patientIdLineEdit)

        # --- Fiducial Table (Collapsible) ---
        self.fiducialTableCollapsible = ctk.ctkCollapsibleButton()
        self.fiducialTableCollapsible.text = "Here a table guiding the position of the needed fiducials nodes for the PHF analysis"
        self.layout.addWidget(self.fiducialTableCollapsible)
        self.fiducialTableLayout = qt.QVBoxLayout(self.fiducialTableCollapsible)

        self.fiducialTable = qt.QTableWidget()
        self.fiducialTable.setRowCount(len(self.fiducialLabels))
        self.fiducialTable.setColumnCount(2)
        self.fiducialTable.setHorizontalHeaderLabels(["ID", "Description"])
        self.fiducialTable.horizontalHeader().setStretchLastSection(True)
        self.fiducialTable.verticalHeader().visible = False
        self.fiducialTable.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)

        for i, (fidID, label) in enumerate(self.fiducialLabels.items()):
            self.fiducialTable.setItem(i, 0, qt.QTableWidgetItem(fidID))
            self.fiducialTable.setItem(i, 1, qt.QTableWidgetItem(label))

        self.fiducialTableLayout.addWidget(self.fiducialTable)


        # --- Preop Section ---
        self.preopGroup = ctk.ctkCollapsibleButton()
        self.preopGroup.text = "Analysis of the Pre-op image"
        self.layout.addWidget(self.preopGroup)
        self.preopLayout = qt.QVBoxLayout(self.preopGroup)
        self.preopGroup.setStyleSheet("background-color: #4CAF50; color: white;")

        self.preopPlaceButton = qt.QPushButton("Place Fiducials (Preop)")
        self.preopPlaceButton.connect('clicked(bool)', lambda: self.onPlaceFiducial("preop"))
        self.preopLayout.addWidget(self.preopPlaceButton)

        self.preopComputeButton = qt.QPushButton("Compute Parameters (Preop)")
        self.preopComputeButton.connect('clicked(bool)', lambda: self.onCalculateMetrics("preop"))
        self.preopLayout.addWidget(self.preopComputeButton)

        self.preopExportButton = qt.QPushButton("Export Fiducials (Preop)")
        self.preopExportButton.connect('clicked(bool)', lambda: self.onExportFiducials("preop"))
        self.preopLayout.addWidget(self.preopExportButton)

        # --- Reset Fiducials Preop Button ---
        self.resetFiducialsButton = qt.QPushButton("Reset Fiducials (Preop)")
        self.resetFiducialsButton.toolTip = "Remove all preop placed fiducial nodes"
        self.resetFiducialsButton.connect('clicked(bool)', self.onResetFiducials_preop)
        self.preopLayout.addWidget(self.resetFiducialsButton)

        # --- Postop Section ---
        self.postopGroup = ctk.ctkCollapsibleButton()
        self.postopGroup.text = "Analysis of the Post-op image"
        self.layout.addWidget(self.postopGroup)
        self.postopLayout = qt.QVBoxLayout(self.postopGroup)
        self.postopGroup.setStyleSheet("background-color: #2196F3; color: white;")

        self.postopPlaceButton = qt.QPushButton("Place Fiducials (Postop)")
        self.postopPlaceButton.connect('clicked(bool)', lambda: self.onPlaceFiducial("postop"))
        self.postopLayout.addWidget(self.postopPlaceButton)

        self.postopComputeButton = qt.QPushButton("Compute Parameters (Postop)")
        self.postopComputeButton.connect('clicked(bool)', lambda: self.onCalculateMetrics("postop"))
        self.postopLayout.addWidget(self.postopComputeButton)

        self.postopExportButton = qt.QPushButton("Export Fiducials (Postop)")
        self.postopExportButton.connect('clicked(bool)', lambda: self.onExportFiducials("postop"))
        self.postopLayout.addWidget(self.postopExportButton)

        # --- Reset Fiducials Preop Button ---
        self.resetFiducialsButton = qt.QPushButton("Reset Fiducials (Postop)")
        self.resetFiducialsButton.toolTip = "Remove all postop placed fiducial nodes"
        self.resetFiducialsButton.connect('clicked(bool)', self.onResetFiducials_postop)
        self.postopLayout.addWidget(self.resetFiducialsButton)

        #Margins of the layout
        self.preopLayout.setContentsMargins(0, 0, 0, 0)
        self.preopLayout.setSpacing(2)
        self.postopLayout.setContentsMargins(0, 0, 0, 0)
        self.postopLayout.setSpacing(2)
        self.layout.setSpacing(3)
        self.layout.setContentsMargins(2, 2, 2, 2)

        self.logic = PHF_analysisLogic()

        # Virtual analysis
        self.vir_analysisComputeButton = qt.QPushButton("Compute virtual humerus reconstruction")
        self.vir_analysisComputeButton.connect('clicked(bool)', lambda: self.onComputeVirtualHumerus("preop"))
        self.layout.addWidget(self.vir_analysisComputeButton)

    """Fonctions definitions"""

    def onInputVolumeChanged(self, volumeNode):
        if volumeNode:
            slicer.util.setSliceViewerLayers(background=volumeNode, fit=True)

    def onLoadDICOM(self):
        slicer.util.selectModule("DICOM")
        qt.QMessageBox.information(slicer.util.mainWindow(), "DICOM Import",
            "The DICOM browser has been opened.\nLoad your X-ray image, then return here to place fiducials.")

    def onResetFiducials_preop(self):
        scene = slicer.mrmlScene

        if not slicer.util.confirmYesNoDisplay("Are you sure you want to remove all preop fiducial nodes?"):
            return

        for node in scene.GetNodesByClass("vtkMRMLMarkupsFiducialNode"):
            name = node.GetName()
            if "preop" in name.lower():
                scene.RemoveNode(node)


    def onResetFiducials_postop(self):
        scene = slicer.mrmlScene

        if not slicer.util.confirmYesNoDisplay("Are you sure you want to remove all preop fiducial nodes?"):
            return

        for node in scene.GetNodesByClass("vtkMRMLMarkupsFiducialNode"):
            name = node.GetName()
            if "postop" in name.lower():
                scene.RemoveNode(node)

    def onPlaceFiducial(self, timepoint):
        """Activate the placement of the fiducial points and label them P1-P13"""
        nodeName = f"Fiducials_{timepoint}"
        self.markupNode = slicer.mrmlScene.GetFirstNodeByName(nodeName)
        if not self.markupNode:
            self.markupNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", nodeName)
            self.markupNode.AddObserver(slicer.vtkMRMLMarkupsNode.PointAddedEvent, self.onFiducialAdded)

        if "preop" in timepoint.lower():
            self.markupNode.GetDisplayNode().SetSelectedColor(0, 1, 0) #green
        elif "postop" in timepoint.lower():
            self.markupNode.GetDisplayNode().SetSelectedColor(0, 0, 1) #blue

        slicer.modules.markups.logic().StartPlaceMode(1)


    def onFiducialAdded(self, caller, event):
        """Automatically label new fiducials based on P1–P13"""
        n = caller.GetNumberOfControlPoints()
        if n <= len(self.fiducialLabels):
            fidName = f"P{n}"
            caller.SetNthControlPointLabel(n - 1, fidName)
        else:
            interactionNode = slicer.app.applicationLogic().GetInteractionNode()
            interactionNode.SetCurrentInteractionMode(slicer.vtkMRMLInteractionNode.ViewTransform)
            slicer.util.showStatusMessage("13 fiducials placed — placement stopped.", 3000)

    def get_input_volume_spacing(self):
        """
        Returns the pixel spacing [spacingX, spacingY, spacingZ]
        from the selected input volume.
        """
        try:
            volumeNode = self.inputVolumeSelector.currentNode()
            if not volumeNode:
                slicer.util.errorDisplay("No volume selected in the input selector.")
                return None

            if not volumeNode.IsA("vtkMRMLScalarVolumeNode"):
                slicer.util.errorDisplay("Selected node is not a scalar volume.")
                return None

            spacing = volumeNode.GetSpacing()  # (x, y, z)
            if spacing == (1.0, 1.0, 1.0):
                slicer.util.warningDisplay("Warning: spacing is default (1.0 mm). Image may lack DICOM geometry.")
            return spacing[:2]  # return only x and y for 2D images

        except Exception as e:
            slicer.util.errorDisplay(f"Error retrieving pixel spacing: {str(e)}")
            return None

    def onExportFiducials(self, timepoint):
        """Export the coordinates of the fiducial points as a .fcsv"""
        nodeName = f"Fiducials_{timepoint}"
        node = slicer.mrmlScene.GetFirstNodeByName(nodeName)
        if not node:
            slicer.util.errorDisplay(f"No fiducials placed for {timepoint}")
            return

        moduleDir = os.path.dirname(os.path.abspath(__file__))
        dataDir = os.path.join(moduleDir, "data")
        os.makedirs(dataDir, exist_ok=True)

        patientID = self.patientIdLineEdit.text.strip()
        if not patientID:
            slicer.util.errorDisplay("Please enter a patient ID before computing.")
            return

        fileName = "fiducials"
        fullFileName = f"{fileName}_{patientID}_{timepoint}.fcsv"

        filePath = os.path.join(dataDir, fullFileName)

        slicer.util.saveNode(node, filePath)

        print(f"Fiducials nodes saved to: {filePath}")


    def onCalculateMetrics(self, timepoint):
        """Computes : GTA - Greater tuberosity angle, NSA - Neck-shaft angle, GTHH - Greater tuberosity to humeral head distance
                    MHL - Medial hinge length, vHH- Virtual humeral head displacement, vGT- Virtual greater tuberosity displacement
                    and save the results as a .csv"""
        nodeName = f"Fiducials_{timepoint}"
        markupNode = slicer.mrmlScene.GetFirstNodeByName(nodeName)
        if not markupNode:
            slicer.util.errorDisplay(f"No fiducials placed for {timepoint}")
            return

        patientID = self.patientIdLineEdit.text.strip()
        if not patientID:
            slicer.util.errorDisplay("Please enter a patient ID.")
            return

        selectedSide  = self.sideSelector.currentText
        pxl_spacing = self.get_input_volume_spacing()
        if not pxl_spacing:
            return

        # Save to temp
        fcsv_path = os.path.join(slicer.app.temporaryPath, f"temp_{timepoint}.fcsv")
        slicer.util.saveNode(markupNode, fcsv_path)
        points = self.logic.loadFCSV(fcsv_path)

        if len(points) < 13:
            slicer.util.errorDisplay(f"Not enough fiducials for {timepoint} (min 13).")
            return

        data = np.array([
            [int(d['label'].lstrip('P')), d['x'], d['y']]
            for d in points
        ])

        # Humeral head circle
        x_values = data[np.isin(data[:, 0], [7, 8, 13, 14, 15]), 1]
        y_values = data[np.isin(data[:, 0], [7, 8, 13, 14, 15]), 2]
        xo, yo, R = self.logic.fitCircle2D(x_values, y_values)
       # self.logic.drawCircle2D((xo, yo), R)
        data_clean = np.delete(data, 0, axis=1)
        # Humeral diaphysis
        humeral_diaphysis = self.logic.compute_humeral_diaphysis(data_clean)
        # Anatomical neck
        NV_points = self.logic.compute_NV_points(data_clean, xo, yo)
        # Intersection of the NSA lines
        intersection_point = self.logic.compute_intersection_point(humeral_diaphysis["DP3"], humeral_diaphysis["DP4"], NV_points["NP3"], NV_points["NP4"])
        # Greater tuberosity tilt
        GT_tilt = self.logic.compute_GT_tilt(humeral_diaphysis["GV1"], humeral_diaphysis["DV2"])
        # Neck shaft angle
        NSA = self.logic.compute_NSA(NV_points["NV2"], humeral_diaphysis["DV1"])
        # Greater tuberosity angle
        GTA = self.logic.compute_GTA(humeral_diaphysis["GV1"], humeral_diaphysis["DV1"])
        # Greater tuberosity - Head distance
        gthh_is, gthh_ml, gthh_3d = self.logic.compute_gthh_displacement(data_clean, selectedSide, pxl_spacing)
        # Medial hinge length
        mhl_is, mhl_ml, mhl_3d = self.logic.compute_mhl_displacement(data_clean, selectedSide, pxl_spacing)

        metrics = {
            "Circle center X": xo,
            "Circle center Y": yo,
            "Radius": R,
            "Intersection point": intersection_point,
            "GT tilt": GT_tilt,
        }

        results = {
            "NSA": NSA,
            "GTA": GTA,
            "GTHH_IS": gthh_is,
            "GTHH_ML": gthh_ml,
            "GTHH_3D": gthh_3d,
            "MHL_IS": mhl_is,
            "MHL_ML": mhl_ml,
            "MHL_3D": mhl_3d
        }
        patientID = f"{self.patientIdLineEdit.text.strip()}_{timepoint}"
        if not patientID:
            slicer.util.errorDisplay("Please enter a patient ID before computing.")
            return

        for key, value in results.items():
            print(f"{key}: {value:.2f}")

        self.logic.saveMetricsToCSV(metrics, "metrics", patientID)
        self.logic.saveMetricsToCSV(humeral_diaphysis, "humeral_diaphysis", patientID)
        self.logic.saveMetricsToCSV(NV_points, "NV_points", patientID)
        self.logic.saveMetricsToCSV(results, "results", patientID)



    def onComputeVirtualHumerus(self, timepoint):
        patientID = self.patientIdLineEdit.text.strip()
        if not patientID:
            slicer.util.errorDisplay("Please enter a patient ID.")
            return

        moduleDir = os.path.dirname(os.path.abspath(__file__))
        dataDir = os.path.join(moduleDir, "data")
        FileName = f"{patientID}_{timepoint}"

        fiducials_name = f"fiducials_{patientID}_preop.fcsv"
        filepath = os.path.join(dataDir, fiducials_name)
        fiducials_dicts = self.logic.loadFCSV(filepath)
        fiducials = np.array([[d['x'], d['y']] for d in fiducials_dicts], dtype=np.float64)

        selectedSide = self.sideSelector.currentText
        pxl_spacing = self.get_input_volume_spacing()
        if not pxl_spacing:
            return
        ratio = 1/pxl_spacing[0]

        Humerus_vir, xo_virt, yo_virt = self.logic.perform_virtual_repositioning(fiducials, FileName, selectedSide)

        translation_cor = fiducials[4] - Humerus_vir[5]

        # Apply translation to points 6–8 and 13–15
        for ipoint in [5, 6, 7, 12, 13, 14]:
            Humerus_vir[ipoint] += translation_cor

        # Translate humeral head center
        xo_virt += translation_cor[0]
        yo_virt += translation_cor[1]

        Humerus_vir, xo_virt, yo_virt = self.logic.perform_virtual_repositioning(fiducials, FileName, selectedSide)
        Humerus_vir = self.logic.reposition_greater_trochanter(Humerus_vir, fiducials, FileName, selectedSide, ratio)


        # Compute direction vectors from results
        humeral_diaphysis = self.logic.load_csv_data("humeral_diaphysis", FileName)
        DV1_str = humeral_diaphysis.loc[humeral_diaphysis['Metric'] == "DV1", 'Value'].values[0]
        DV1 = np.array(ast.literal_eval(DV1_str), dtype=np.float64)

        DV2_str = humeral_diaphysis.loc[humeral_diaphysis['Metric'] == "DV2", 'Value'].values[0]
        DV2 = np.array(ast.literal_eval(DV2_str), dtype=np.float64)

        NV_points = self.logic.load_csv_data("NV_points", FileName)
        NV2_str = NV_points.loc[NV_points['Metric'] == "NV2", 'Value'].values[0]
        NV2 = np.array(ast.literal_eval(NV2_str), dtype=np.float64)

        # Set up
        side = selectedSide
        H_vir = Humerus_vir
        H_pre = fiducials

        # GT superior medial line
        GP1 = H_vir[12]
        GP2 = H_vir[9]
        GV1 = (GP2 - GP1) / np.linalg.norm(GP2 - GP1)
        GP3 = 200 * GV1 + GP1
        GP4 = -250 * GV1 + GP1

        # GT Tilt
        dot_GT = np.dot(GV1, DV2) / (np.linalg.norm(GV1) * np.linalg.norm(DV2))
        GT_angle = np.degrees(np.arccos(np.clip(dot_GT, -1.0, 1.0)))
        GT = 180 - GT_angle if GT_angle > 90 else GT_angle

        # NSA
        dot_NSA = np.dot(NV2, DV1) / (np.linalg.norm(NV2) * np.linalg.norm(DV1))
        NSA_angle = np.degrees(np.arccos(np.clip(dot_NSA, -1.0, 1.0)))
        NSA = 180 - NSA_angle if NSA_angle < 90 else NSA_angle

        # GTA
        dot_GTA = np.dot(GV1, DV1) / (np.linalg.norm(GV1) * np.linalg.norm(DV1))
        GTA = np.degrees(np.arccos(np.clip(dot_GTA, -1.0, 1.0)))

        # GT–HH displacement
        GTHH_IS = -(H_vir[9, 1] - H_vir[8, 1]) / ratio
        GTHH_ML = (H_vir[9, 0] - H_vir[8, 0]) / ratio
        if side == "Right":
            GTHH_ML = -GTHH_ML
        GTHH_3D = np.linalg.norm(H_vir[9] - H_vir[8]) / ratio

        # Medial hinge length
        MHL_IS = -(H_vir[6, 1] - H_vir[5, 1]) / ratio
        MHL_ML = (H_vir[6, 0] - H_vir[5, 0]) / ratio
        if side == "Right":
            MHL_ML = -MHL_ML
        MHL_3D = np.linalg.norm(H_vir[6] - H_vir[5]) / ratio

        # Virtual GT displacement
        vGT_IS = -(H_pre[9, 1] - H_vir[8, 1]) / ratio
        vGT_ML = (H_pre[9, 0] - H_vir[8, 0]) / ratio
        if side == "Right":
            vGT_ML = -vGT_ML
        vGT_3D = np.linalg.norm(H_pre[9] - H_vir[8]) / ratio

        # Virtual HH displacement
        vHH_IS = -(H_pre[8, 1] - H_vir[8, 1]) / ratio
        vHH_ML = (H_pre[8, 0] - H_vir[8, 0]) / ratio
        if side == "Right":
            vHH_ML = -vHH_ML
        vHH_3D = np.linalg.norm(H_pre[8] - H_vir[8]) / ratio

        # Store everything
        outcomes = {
            "GT": GT,
            "NSA": NSA,
            "GTA": GTA,
            "GTHH_IS": GTHH_IS,
            "GTHH_ML": GTHH_ML,
            "GTHH_3D": GTHH_3D,
            "MHL_IS": MHL_IS,
            "MHL_ML": MHL_ML,
            "MHL_3D": MHL_3D,
            "vGT_IS": vGT_IS,
            "vGT_ML": vGT_ML,
            "vGT_3D": vGT_3D,
            "vHH_IS": vHH_IS,
            "vHH_ML": vHH_ML,
            "vHH_3D": vHH_3D
        }

        for key, value in outcomes.items():
            print(f"{key}: {value:.2f}")

        output_file = os.path.join(dataDir, f"outcomes_{FileName}.csv")
        with open(output_file, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Metric", "Value"])
            for key, value in outcomes.items():
                writer.writerow([key, value])

#
# PHF_analysisLogic
#


class PHF_analysisLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def getFiducialPoints(self, markupNode):
        coords = []
        for i in range(markupNode.GetNumberOfFiducials()):
            pos = [0, 0, 0]
            markupNode.GetNthFiducialPosition(i, pos)
            coords.append(pos[:2])
        return np.array(coords)

    def loadFCSV(self, filePath):
        points = []
        with open(filePath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split(',')
                if len(parts) < 12:
                    continue
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                label = parts[11]
                points.append({'label': label, 'x': x, 'y': y, 'z': z})
        return points

    def fitCircle2D(self, x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        A = np.c_[2 * x, 2 * y, np.ones(len(x))]
        b = x ** 2 + y ** 2
        coeffs = np.linalg.lstsq(A, b, rcond=None)[0]
        xo, yo = coeffs[0], coeffs[1]
        R = np.sqrt(coeffs[2] + xo ** 2 + yo ** 2)
        return xo, yo, R

    def compute_humeral_diaphysis(self, AL):
        """
        Compute humeral diaphysis points and vectors from anatomical landmarks.
        """

        # ---- Diaphyseal points ----
        DP1 = (AL[0] + AL[3]) / 2
        DP2 = (AL[1] + AL[2]) / 2
        DV1 = (DP2 - DP1)
        DV1 = DV1 / np.linalg.norm(DV1)  # Tangent vector
        DV2 = np.array([DV1[1], -DV1[0]])  # Normal vector in XY plane

        DP3 = DP1 + 600 * DV1
        DP4 = DP1 - 500 * DV1
        DP5 = AL[11] + 100 * DV2
        DP6 = AL[11] - 250 * DV2

        # ---- Great Tuberosity line ----
        GP1 = AL[11]
        GP2 = AL[8]
        GV1 = (GP2 - GP1)
        GV1 = GV1 / np.linalg.norm(GV1)

        GP3 = GP1 + 200 * GV1
        GP4 = GP1 - 250 * GV1

        return {
            "DP1": DP1, "DP2": DP2,
            "DP3": DP3, "DP4": DP4,
            "DP5": DP5, "DP6": DP6,
            "DV1": DV1, "DV2": DV2,
            "GP1": GP1, "GP2": GP2,
            "GP3": GP3, "GP4": GP4,
            "GV1": GV1
        }

    def compute_NV_points(self, AL, xo, yo):
        """
        Calcule les points et vecteurs normaux basés sur les points 7 et 8 (NP1, NP2)
        ainsi que leur base tangente et normale projetée depuis le centre du cercle.
        """
        NP1 = AL[6, :2]
        NP2 = AL[7, :2]

        NV1 = NP2 - NP1
        NV1 = NV1 / np.linalg.norm(NV1)
        NV2 = np.array([NV1[1], -NV1[0]])

        center = np.array([xo, yo])

        NP3 = 150 * NV2 + center
        NP4 = -500 * NV2 + center

        return {
            "NP1": NP1,
            "NP2": NP2,
            "NV1": NV1,
            "NV2": NV2,
            "NP3": NP3,
            "NP4": NP4
        }

    def compute_intersection_point(self, DP3, DP4, NP3, NP4):
        """
        Calcule le point d'intersection moyen des deux droites formées par DP3-DP4 et NP3-NP4.
        """
        x1, y1 = DP3[0], DP3[1]
        x2, y2 = np.ravel(DP4)
        x3, y3 = np.ravel(NP3)
        x4, y4 = np.ravel(NP4)

        denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

        if denom == 0:
            raise ValueError("Les droites sont parallèles, pas d'intersection.")

        u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / denom
        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom

        xi = ((x3 + u * (x4 - x3)) + (x1 + t * (x2 - x1))) / 2
        yi = ((y3 + u * (y4 - y3)) + (y1 + t * (y2 - y1))) / 2

        return np.array([xi, yi])

    def compute_GT_tilt(self, GV1, DV2):
        # Assure-toi que ce sont bien des vecteurs 1D
        GV1 = np.ravel(GV1)
        DV2 = np.ravel(DV2)

        dot_product = np.dot(GV1, DV2)
        norm_product = np.linalg.norm(GV1) * np.linalg.norm(DV2)

        angle_rad = math.acos(np.clip(dot_product / norm_product, -1.0, 1.0))
        angle_deg = math.degrees(angle_rad)

        if angle_deg > 90:
            angle_deg = 180 - angle_deg

        return angle_deg

    def compute_NSA(self, NV2, DV1):
        NV2 = np.ravel(NV2)
        DV1 = np.ravel(DV1)

        dot_product = np.dot(NV2, DV1)
        norm_product = np.linalg.norm(NV2) * np.linalg.norm(DV1)

        angle_rad = math.acos(np.clip(dot_product / norm_product, -1.0, 1.0))
        angle_deg = math.degrees(angle_rad)

        if angle_deg < 90:
            angle_deg = 180 - angle_deg

        return angle_deg

    def compute_GTA(self, GV1, DV1):
        # Greater tuberosity angle (GTA)
        dot_product = np.dot(GV1, DV1)
        norm_product = np.linalg.norm(GV1) * np.linalg.norm(DV1)
        GTA = math.degrees(math.acos(np.clip(dot_product / norm_product, -1.0, 1.0)))

        return GTA

    def compute_gthh_displacement(self, AL, side, pxl_spacing):
        """
        Compute GTHH inferosuperior, mediolateral, and 3D displacements.
        """
        AL = np.asarray(AL)
        ratio = 1/pxl_spacing[0]

        gthh_is = -(AL[8][1] - AL[7][1]) / ratio

        delta_ml = AL[8][0] - AL[7][0]
        if side == "Right":
            gthh_ml = -delta_ml / ratio
        elif side == "Left":
            gthh_ml = delta_ml / ratio
        else:
            raise ValueError("Invalid side. Must be 'Right' or 'Left'.")

        delta_x = AL[8][0] - AL[7][0]
        delta_y = AL[8][1] - AL[7][1]
        gthh_3d = np.sqrt(delta_x ** 2 + delta_y ** 2) / ratio

        return gthh_is, gthh_ml, gthh_3d

    def compute_mhl_displacement(self, AL, side, pxl_spacing):
        """
        Compute MHL (Medial Hinge Length) displacements:
        - Inferosuperior (IS)
        - Mediolateral (ML)
        - 3D (combined)
        """
        ratio = 1 / pxl_spacing[0]

        delta_y = AL[5][1] - AL[4][1]
        delta_x = AL[5][0] - AL[4][0]

        mhl_is = -delta_y / ratio

        if side == "Right":
            mhl_ml = -delta_x / ratio
        elif side == "Left":
            mhl_ml = delta_x / ratio
        else:
            raise ValueError("Invalid side: must be 'Right' or 'Left'")

        mhl_3d = np.sqrt(delta_x ** 2 + delta_y ** 2) / 7  # Hardcoded 7 in original MATLAB

        return mhl_is, mhl_ml, mhl_3d

    def saveMetricsToCSV(self, metrics, fileName, patientID):
        # Get the path to this script's directory (i.e., the module directory)
        moduleDir = os.path.dirname(os.path.abspath(__file__))
        dataDir = os.path.join(moduleDir, "data")
        os.makedirs(dataDir, exist_ok=True)

        fullFileName = f"{fileName}_{patientID}.csv"

        filePath = os.path.join(dataDir, fullFileName)

        with open(filePath, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Metric", "Value"])
            for key, value in metrics.items():
                if isinstance(value, (list, tuple)):
                    value_str = ",".join(str(v) for v in value)
                elif hasattr(value, "shape"):  # likely a numpy array
                    value_str = ",".join(str(v) for v in value.flatten())
                else:
                    value_str = str(value)
                writer.writerow([key, value_str])

        print(f"File saved as: {filePath}")

    def load_csv_data(self, metric, id_op):

        moduleDir = os.path.dirname(os.path.abspath(__file__))
        dataDir = os.path.join(moduleDir, "data")

        FileName = f"{metric}_{id_op}.csv"
        filePath = os.path.join(dataDir, FileName)
        data = pd.read_csv(filePath)

        return data

    # ----- Virtual Humerus -----

    def perform_virtual_repositioning(self, fiducials, FileName, selectedSide):

        # Constants
        NSA_ref = float(135)

        results = self.load_csv_data("results", FileName)
        metrics = self.load_csv_data("metrics", FileName)

        # Extract required values
        xo = float(metrics.loc[metrics['Metric'] == "Circle center X", 'Value'].values[0])
        yo = float(metrics.loc[metrics['Metric'] == "Circle center Y", 'Value'].values[0])
        R = float(metrics.loc[metrics['Metric'] == "Radius", 'Value'].values[0])
        pre_NSA = float(results.loc[results['Metric'] == "NSA", 'Value'].values[0])

        AL_pre = fiducials
        AL_virt = np.full((15, 2), np.nan)

        # Copy unaffected points
        AL_virt[[0, 1, 2, 3, 4, 10], :] = AL_pre[[0, 1, 2, 3, 4, 10], :]

        # Determine angle correction based on side
        if selectedSide == "Right":
            angle_cor = pre_NSA - NSA_ref
        elif selectedSide == "Left":
            angle_cor = -(pre_NSA - NSA_ref)
        else:
            raise ValueError("Invalid side: must be 'Right' or 'Left'")

        # Rotation matrix for angular correction
        angle_rad = np.radians(angle_cor)
        rotation_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
        ])

        for ipoint in [5, 6, 7, 12]: # /!\ 13, 14 ???? Il n'y que 13 points...
            point = AL_pre[ipoint] - np.array([xo, yo])
            AL_virt[ipoint] = (rotation_matrix @ point) + np.array([xo, yo])

        temp_center = rotation_matrix @ np.array([xo, yo])
        xo_virt = temp_center[0]
        yo_virt = temp_center[1]

        translation_cor = AL_pre[4] - AL_virt[5]

        for ipoint in [5, 6, 7, 12, 13, 14]:
            AL_virt[ipoint] += translation_cor

        xo_virt += translation_cor[0]
        yo_virt += translation_cor[1]

        return AL_virt, xo_virt, yo_virt


    def reposition_greater_trochanter(self, Humerus_vir, Humerus_pre, FileName, selectedSide, ratio):
        GTA_ref = 65.2
        results = self.load_csv_data("results", FileName)
        pre_GTA = float(results.loc[results['Metric'] == "GTA", 'Value'].values[0])

        # Step 1: Determine correction angle
        if selectedSide == "Right":
            angle_cor = -(pre_GTA - GTA_ref)
        elif selectedSide == "Left":
            angle_cor = pre_GTA - GTA_ref
        else:
            raise ValueError("Invalid side: must be 'Right' or 'Left'")

        # Rotation matrix around point 9 (index 8)
        center = Humerus_pre[8]
        rotation_matrix = np.array([
            [np.cos(np.radians(angle_cor)), -np.sin(np.radians(angle_cor))],
            [np.sin(np.radians(angle_cor)), np.cos(np.radians(angle_cor))]
        ])

        # Step 1: Rotate points 9,10,12 (indices 8,9,11)
        for ipoint in [8, 9, 11]:
            rel = Humerus_pre[ipoint] - center
            Humerus_vir[ipoint] = (rotation_matrix @ rel) + center

        # Step 2: Translate to match point 8 with pre position
        translation_cor = Humerus_vir[7] - Humerus_vir[8]  # point 8 - point 9

        for ipoint in [8, 9, 11]:
            Humerus_vir[ipoint] += translation_cor

        # Step 3: If initial offset between AL(8) and AL(9) < 1mm, reapply simple rotation
        dist = np.linalg.norm(Humerus_pre[7] - Humerus_pre[8]) / ratio
        if dist < 1:
            # Re-rotate around (xo, yo) instead
            xo, yo = Humerus_pre[8]
            center = np.array([xo, yo])
            for ipoint in [8, 9, 11]:
                rel = Humerus_pre[ipoint] - center
                Humerus_vir[ipoint] = (rotation_matrix @ rel) + center

            # Reapply translation again
            translation_cor = Humerus_vir[7] - Humerus_vir[8]
            for ipoint in [8, 9, 11]:
                Humerus_vir[ipoint] += translation_cor
        else:
            # Use numerical optimization (as in MATLAB's fminsearch)
            def error_fn(x):
                angle, dx, dy = x
                R = np.array([
                    [np.cos(np.radians(angle)), -np.sin(np.radians(angle))],
                    [np.sin(np.radians(angle)), np.cos(np.radians(angle))]
                ])
                d = np.array([dx, dy])
                p1 = (R @ Humerus_pre[8]) + d
                p2 = (R @ Humerus_pre[9]) + d
                err = np.linalg.norm(p1 - Humerus_vir[7]) + np.linalg.norm(p2 - Humerus_vir[10])
                return err

            x0 = [2, 1, 1]
            res = minimize(error_fn, x0, method='Nelder-Mead')

            angle_opt, dx_opt, dy_opt = res.x
            R = np.array([
                [np.cos(np.radians(angle_opt)), -np.sin(np.radians(angle_opt))],
                [np.sin(np.radians(angle_opt)), np.cos(np.radians(angle_opt))]
            ])
            d = np.array([dx_opt, dy_opt])

            # Debug print for mm errors (optional)
            ecart_prox = np.linalg.norm((R @ Humerus_pre[8] + d) - Humerus_vir[7]) / ratio
            ecart_dist = np.linalg.norm((R @ Humerus_pre[9] + d) - Humerus_vir[10]) / ratio
            print(f"proximal drift: {ecart_prox:.2f} mm")
            print(f"distal drift: {ecart_dist:.2f} mm")

            # Apply optimized transform
            for ipoint in [8, 9, 11]:
                Humerus_vir[ipoint] = (R @ Humerus_pre[ipoint]) + d

        return Humerus_vir


#
# PHF_analysisTest
#


class PHF_analysisTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """Do whatever is needed to reset the state - typically a scene clear will be enough."""
        slicer.mrmlScene.Clear()

    def runTest(self):
        """Run as few or as many tests as needed here."""
        self.setUp()
        self.test_PHF_analysis1()

    def test_PHF_analysis1(self):
        """Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("Starting the test")

        # Get/create input data

        import SampleData

        registerSampleData()
        inputVolume = SampleData.downloadSample("PHF_analysis1")
        self.delayDisplay("Loaded test data set")

        inputScalarRange = inputVolume.GetImageData().GetScalarRange()
        self.assertEqual(inputScalarRange[0], 0)
        self.assertEqual(inputScalarRange[1], 695)

        outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
        threshold = 100

        # Test the module logic

        logic = PHF_analysisLogic()

        # Test algorithm with non-inverted threshold
        logic.process(inputVolume, outputVolume, threshold, True)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], threshold)

        # Test algorithm with inverted threshold
        logic.process(inputVolume, outputVolume, threshold, False)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], inputScalarRange[1])

        self.delayDisplay("Test passed")
