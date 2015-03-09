"""
    Skeletonizer: Python Cell Morphology Analysis and Construction Toolkit

    KAUST, BESE, Neuro-Inspired Computing Project
    (c) 2014-2015. All rights reserved.
"""

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'skeletonizer',
    'description': 'Cell morphology analysis and construction tool for BBPSDK',
    'long_description': ''' 
                        Program for constructing a BBPSDK cell morphology (into HDF5 format) from
                        mesh (Blender source, exported into VRML for import into Avizo) and 
                        skeletonization data (Avizo Amiramesh ASCII format).
                        Generates accurate cross-sectional data from mesh and skeleton points (CSV format).
                        '''
    'author': 'Neuro-Inspired Computing Team: Glendon Holst, Heikki Lehvaslaiho, et. al.',
    'author_email': 'glendon.holst@kaust.edu.sa',
    'maintainer': 'Neuro-Inspired Computing Team: Daniya Boges'
    'maintainer_email': 'daniya.boges@kaust.edu.sa'
    'url': 'https://bitbucket.org/holstgr/skeletonizer',
    'version': '1.0.0b1',
    'license': 'MIT',
    'install_requires': ['unittest','bbp'],
    'packages': ['skeletonizer'],
    'py_module': ['skeletonizer.amiramesh',
                  'skeletonizer.graphs',
                  'skeletonizer.maths',
                  'skeletonizer.morphology'
                 ]
    'scripts': ['bin/skeletonize.py', 'bin/skeleton_annotate.py'],
    'data_files': [('test',['data/test.blend',
                            'data/test.SptGraph.am',
                            'data/test.SptGraph.annotations.json'
                           ]
                  )],
    'keywords': ('morphology', 'analysis', 'generation', 'skeletonization'),
    'classifiers': ['License :: OSI Approved :: MIT License',
                    'Development Status :: 4 - Beta',
                    'Programming Language :: Python :: 2.7',
                    'Natural Language :: English',
                    'Environment :: Console',
                    'Operating System :: POSIX :: Linux',
                    'Topic :: Scientific/Engineering :: Bio-Informatics',
                    'Intended Audience :: Science/Research'
                   ]

}
setup(**config)
