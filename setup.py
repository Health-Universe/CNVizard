from setuptools import setup, find_packages

setup(
    name='CNVizard',
    version='0.1',
    description='A streamlit app for visualizing germline copy number variants',
    author='Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft',
    author_email='jerkrause@ukaachen.de',
    url='https://github.com/jerkrause/CNVizard',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'streamlit',
        'pandas',
        'numpy',
        'xlsxwriter',
        'pyarrow',
        'plotly',
        'cnvlib',
        'python-dotenv'
    ],
    entry_points={
        'console_scripts': [
            'run_cnvizard=run_visualizer:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
