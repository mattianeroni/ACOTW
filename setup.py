from setuptools import setup

setup(
    name='ACOTW',
    version='0.0.1',
    description='Ant Colony Optimization for path finding avoiding collisions between AGVs.',
    author='Mattia Neroni, Ph.D, Eng.',
    author_email='mattia.neroni@ahead-research.com',
    url='https://github.com/mattianeroni/ACOTW',
    packages=[
        "acotw"
    ],
    python_requires='>=3.9',
    classifiers=[
        "Development Status :: 3 - Alpha"
    ]
)