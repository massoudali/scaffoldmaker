"""
Generates a solid Bone using a ShieldMesh of all cube elements,
 with variable numbers of elements in major, minor, shell and axial directions.
"""

from __future__ import division

import copy

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.field import Field
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.cylindermesh import CylinderMesh, CylinderShape, CylinderEnds, CylinderCentralPath
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from scaffoldmaker.utils.derivativemoothing import DerivativeSmoothing
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm


class MeshType_3d_solidbone2 (Scaffold_base):
    """
Generates a solid cylinder using a ShieldMesh of all cube elements,
with variable numbers of elements in major, minor, shell and axial directions.
    """
    centralPathDefaultScaffoldPackages = {
        'Cylinder 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 3.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 2.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
                ])
        })
    }

    @staticmethod
    def getName():
        return '3D Bone'

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        centralPathOption = cls.centralPathDefaultScaffoldPackages['Cylinder 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            'Number of elements across shell': 0,
            'Number of elements across transition': 1,
            'Number of elements along': 4,
            'Shell element thickness proportion': 1.0,
            'Lower half': False,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements across major': 1,
            'Refine number of elements along': 1,
            'Upper scale': 1.0,
            'Upper scale_Z': 1.0,
            'Lower scale': 1.0,
            'Lower scale_Z': 1.0
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements across major',
            'Number of elements across minor',
            # 'Number of elements across shell',
            # 'Number of elements across transition',
            'Number of elements along',
            # 'Shell element thickness proportion',
            # 'Lower half',
            'Refine',
            'Refine number of elements across major',
            'Refine number of elements along',
            'Upper scale',
            'Upper scale_Z',
            'Lower scale',
            'Lower scale_Z'
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        dependentChanges = False

        if options['Number of elements across major'] < 4:
            options['Number of elements across major'] = 4
        if options['Number of elements across major'] % 2:
            options['Number of elements across major'] += 1

        if options['Number of elements across minor'] < 4:
            options['Number of elements across minor'] = 4
        if options['Number of elements across minor'] % 2:
            options['Number of elements across minor'] += 1
        if options['Number of elements along'] < 1:
            options['Number of elements along'] = 1
        if options['Number of elements across transition'] < 1:
            options['Number of elements across transition'] = 1
        Rcrit = min(options['Number of elements across major']-4, options['Number of elements across minor']-4)//2
        if options['Number of elements across shell'] + options['Number of elements across transition'] - 1 > Rcrit:
            dependentChanges = True
            options['Number of elements across shell'] = Rcrit
            options['Number of elements across transition'] = 1

        if options['Shell element thickness proportion'] < 0.15:
            options['Shell element thickness proportion'] = 1.0

        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        centralPath = options['Central path']
        full = not options['Lower half']
        elementsCountAcrossMajor = options['Number of elements across major']
        if not full:
            elementsCountAcrossMajor //= 2
        elementsCountAcrossMinor = options['Number of elements across minor']
        elementsCountAcrossShell = options['Number of elements across shell']
        elementsCountAcrossTransition = options['Number of elements across transition']
        elementsCountAlong = options['Number of elements along']
        shellProportion = options['Shell element thickness proportion']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        cylinderCentralPath = CylinderCentralPath(region, centralPath, elementsCountAlong)

        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL if full else CylinderShape.CYLINDER_SHAPE_LOWER_HALF

        base = CylinderEnds(elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                            elementsCountAcrossTransition,
                            shellProportion,
                            [0.0, 0.0, 0.0], cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
                            cylinderCentralPath.minorRadii[0])
        cylinder1 = CylinderMesh(fm, coordinates, elementsCountAlong, base,
                                 cylinderShape=cylinderShape,
                                 cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fm.findMeshByDimension(3)

        sphere_shape = SphereShape.SPHERE_SHAPE_FULL
        sphere_base = cylinder1._ellipses[0]
        sphere_centre = sphere_base.centre
        sphere_radius_3 = options['Lower scale_Z']
        axes = [sphere_base.majorAxis, sphere_base.minorAxis,
                vector.setMagnitude(vector.crossproduct3(sphere_base.majorAxis, sphere_base.minorAxis),
                                    sphere_radius_3)]
        elementsCountAcross = [cylinder1._elementsCountAcrossMajor, cylinder1._elementsCountAcrossMinor, 4]
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [0, -1]]
        sphereBoxDerivatives = [1, 3, 2]

        sphere1 = SphereMesh(fm, coordinates, [0.0, 0.0, 0.0], axes, elementsCountAcross,
                             0, 1, 1.0,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False)

        hemisphere = ShieldMesh3D(elementsCountAcross, 0)

        # get hemisphere nodes from both cylinder end and top of the sphere and mix them
        btx = hemisphere.px
        btd1 = hemisphere.pd1
        btd2 = hemisphere.pd2
        btd3 = hemisphere.pd3

        hemisphere._boxDerivatives = sphere1._shield3D._boxDerivatives
        hemisphere._boxMapping = sphere1._shield3D._boxMapping
        hemisphere._box_deriv_mapping = sphere1._shield3D._box_deriv_mapping
        hemisphere._element_needs_scale_factor = sphere1._shield3D._element_needs_scale_factor
        hemisphere._xi_mapping = sphere1._shield3D._xi_mapping
        hemisphere._xi_signs = sphere1._shield3D._xi_signs

        for n3 in range(elementsCountAcross[2] + 1):
            for n2 in range(elementsCountAcross[0] + 1):
                for n1 in range(elementsCountAcross[1] + 1):
                    if n3 < elementsCountAcross[2] // 2:
                        if sphere1._shield3D.px[n3][n2][n1]:
                            # hemisphere.nodeId[n3][n2][n1] = sphere1._shield3D.nodeId[n3][n2][n1]
                            btx[n3][n2][n1] = vector.addVectors([sphere1._shield3D.px[n3][n2][n1], sphere_centre],
                                                                [1, 1])
                            btd1[n3][n2][n1] = sphere1._shield3D.pd1[n3][n2][n1]
                            btd2[n3][n2][n1] = sphere1._shield3D.pd2[n3][n2][n1]
                            btd3[n3][n2][n1] = sphere1._shield3D.pd3[n3][n2][n1]

                    # cylinder end
                    elif n3 == elementsCountAcross[2] // 2:
                        # find nodes on the triple line. Note that cylinder and sphere have a little bit different
                        # numbering for nodes on the triple line
                        n2c, n1c = n2, n1
                        if n2 < 1 and n1 == n2:
                            n1c = 1
                        elif n2 < 1 and n1 == elementsCountAcross[1] - n2:
                            n1c = elementsCountAcross[1] - 1
                        elif n2 > elementsCountAcross[1] - 1:
                            if n1 == elementsCountAcross[1] - n2:
                                n1c = 1
                            elif n1 == n2:
                                n1c = elementsCountAcross[1] - 1
                        hemisphere.nodeId[n3][n2][n1] = cylinder1._shield.nodeId[0][n2c][n1c]
#******************************************************************************************************************************
        # generate hemisphere extra nodes.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [0, 1]]
        radius_scale = options['Lower scale']
        for n2 in range(5):
            for n1 in range(5):
                if hemisphere.px[0][n2][n1]:
                    x = hemisphere.px[0][n2][n1]
                    hemisphere.px[0][n2][n1] = [radius_scale * x[0], radius_scale * x[1], x[2]]
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        nodeIdentifier = hemisphere.generateNodes(fm, coordinates, nodeIdentifier,
                                                  rangeOfRequiredElements)
        # rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [0, 1]]
        # nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        # nodeIdentifier = hemisphere.generateNodes(fm, coordinates, nodeIdentifier,
        #                                           rangeOfRequiredElements)
#**********************************************************************************************************************


        # generate hemisphere elements.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [0, 2]]
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        elementIdentifier = hemisphere.generateElements(fm, coordinates, elementIdentifier,
                                                        rangeOfRequiredElements)

        sphere_base = cylinder1._ellipses[-1]
        sphere_centre = sphere_base.centre
        sphere_radius_3 = options['Upper scale_Z']
        axes = [sphere_base.majorAxis, sphere_base.minorAxis,
                vector.setMagnitude(vector.crossproduct3(sphere_base.majorAxis, sphere_base.minorAxis),
                                    sphere_radius_3)]
        elementsCountAcross = [cylinder1._elementsCountAcrossMajor, cylinder1._elementsCountAcrossMinor, 4]
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [0, -1]]
        sphereBoxDerivatives = [1, 3, 2]

        sphere2 = SphereMesh(fm, coordinates, [0.0, 0.0, 0.0], axes, elementsCountAcross,
                             0, 1, 1.0,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False)

        hemisphere = ShieldMesh3D(elementsCountAcross, 0)

        # get hemisphere nodes from both cylinder end and top of the sphere and mix them
        btx = hemisphere.px
        btd1 = hemisphere.pd1
        btd2 = hemisphere.pd2
        btd3 = hemisphere.pd3

        hemisphere._boxDerivatives = sphere2._shield3D._boxDerivatives
        hemisphere._boxMapping = sphere2._shield3D._boxMapping
        hemisphere._box_deriv_mapping = sphere2._shield3D._box_deriv_mapping
        hemisphere._element_needs_scale_factor = sphere2._shield3D._element_needs_scale_factor
        hemisphere._xi_mapping = sphere2._shield3D._xi_mapping
        hemisphere._xi_signs = sphere2._shield3D._xi_signs

        for n3 in range(elementsCountAcross[2] + 1):
            for n2 in range(elementsCountAcross[0] + 1):
                for n1 in range(elementsCountAcross[1] + 1):
                    if n3 > elementsCountAcross[2] // 2:
                        if sphere2._shield3D.px[n3][n2][n1]:
                            # hemisphere.nodeId[n3][n2][n1] = sphere1._shield3D.nodeId[n3][n2][n1]
                            btx[n3][n2][n1] = vector.addVectors([sphere2._shield3D.px[n3][n2][n1], sphere_centre],
                                                                [1, 1])
                            btd1[n3][n2][n1] = sphere2._shield3D.pd1[n3][n2][n1]
                            btd2[n3][n2][n1] = sphere2._shield3D.pd2[n3][n2][n1]
                            btd3[n3][n2][n1] = sphere2._shield3D.pd3[n3][n2][n1]

                    # cylinder end
                    elif n3 == elementsCountAcross[2] // 2:
                        # find nodes on the triple line. Note that cylinder and sphere have a little bit different
                        # numbering for nodes on the triple line
                        n2c, n1c = n2, n1
                        if n2 < 1 and n1 == n2:
                            n1c = 1
                        elif n2 < 1 and n1 == elementsCountAcross[1] - n2:
                            n1c = elementsCountAcross[1] - 1
                        elif n2 > elementsCountAcross[1] - 1:
                            if n1 == elementsCountAcross[1] - n2:
                                n1c = 1
                            elif n1 == n2:
                                n1c = elementsCountAcross[1] - 1
                        hemisphere.nodeId[n3][n2][n1] = cylinder1._shield.nodeId[-1][n2c][n1c]
#******************************************************************************************************************
        # generate hemisphere extra nodes.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [3, 4]]
        radius_scale = options['Upper scale']
        for n2 in range(5):
            for n1 in range(5):
                if hemisphere.px[4][n2][n1]:
                    x = hemisphere.px[4][n2][n1]
                    hemisphere.px[4][n2][n1] = [radius_scale * x[0], radius_scale * x[1], x[2]]
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        nodeIdentifier = hemisphere.generateNodes(fm, coordinates, nodeIdentifier,
                                                  rangeOfRequiredElements)
        # rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [3, 4]]
        # nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        # nodeIdentifier = hemisphere.generateNodes(fm, coordinates, nodeIdentifier,
        #                                           rangeOfRequiredElements)

#************************************************************************************************************************
        # generate hemisphere elements.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [2, 4]]
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        elementIdentifier = hemisphere.generateElements(fm, coordinates, elementIdentifier,
                                                        rangeOfRequiredElements)

        # AnnotationGroup.createMarkerNode(startNodeIdentifier=1, materialCoordinatesField: FieldFiniteElement = None, materialCoordinates = None, element = None, xi = [0.0, 0.0, 0.0])


        # markerTermNameBoneCoordinatesMap = {
        #     'tibial tuberosity': [-5.076472492200136e+01, -4.592226612078402e+01, -1.261953033704384e+03],
        # }
        # annotationGroups = []
        # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        # print(nodeIdentifier)
        # for termName, boneCoordinatesValues in markerTermNameBoneCoordinatesMap.items():
        #     annotationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, ('tibia point', 'FM:63'))
        #     annotationGroup.createMarkerNode(nodeIdentifier, coordinates, boneCoordinatesValues)
        #     nodeIdentifier += 1



        smoothing = DerivativeSmoothing(region, coordinates)
        smoothing.smooth(True)

        annotationGroup = []
        return annotationGroup

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCountAcrossMajor = options['Refine number of elements across major']
        refineElementsCountAlong = options['Refine number of elements along']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCountAcrossMajor, refineElementsCountAlong, refineElementsCountAcrossMajor)

