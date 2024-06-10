from setuptools import setup, find_packages

setup(
    name="CNVizard",
    version="0.1",
    description="A tool for visualizing germline copy number variants",
    author="Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft",
    author_email="jerkrause@ukaachen.de",
    url="https://github.com/IHGGM-Aachen/CNVizard",
    packages=find_packages(include=['cnvizard', 'cnvizard.*']),
    install_requires=[
        "streamlit",
        "pandas",
        "numpy",
        "pyarrow",
        "CNVkit",
        "fastparquet",
        "plotly",
        "python-dotenv"
    ],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "cnvizard=cnvizard.app:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
