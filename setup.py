#!/usr/bin/env python

from distutils.core import setup

def main():
    setup(name='earm2',
          version=1.0,
          description='Extrinsic Apoptosis Reaction Model 2.0',
          author='?',
          author_email='?',
          packages=['earm2'],
          requires=['pysb'],
          keywords=['systems', 'biology', 'model', 'rules'],
          classifiers=[
            'Intended Audience :: Science/Research',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            ],
          )

if __name__ == '__main__':
    main()
