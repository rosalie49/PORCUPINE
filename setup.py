from setuptools import setup, find_packages

setup(
    name='PORCUPINE',
    version='1.0.0',
    description='A package implementing a Principal Components Analysis (PCA)-based approach to identify biological pathways driving inter-tumour heterogeneity in gene regulatory networks.',
    url='https://github.com/rosalie49/PORCUPINE',
    author='Rosalie Coulon',
    author_email='rosalie.coulon@hotmail.fr',
    maintainer='Ladislav Hovan',
    maintainer_email='ladislav.hovan@ncmm.uio.no',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'pyreadr',
        'pandas',
        'numpy',
        'scikit-learn',
        'scipy',
        'yellowbrick',
        'seaborn',
        'matplotlib',
    ],
    zip_safe=False
)