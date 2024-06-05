from setuptools import setup, find_packages

setup(
    name="CNVizard",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "streamlit",
        "pandas",
        "numpy",
        "xlsxwriter",
        "pyarrow",
        "plotly",
        "cnvkit" 
    ],
    entry_points={
        'console_scripts': [
            'run_visualizer=cnvizard.run_visualizer:main',
            'create_reference_files=cnvizard.reference_builder.create_reference_files:main',
            'merge_reference_files=cnvizard.reference_builder.merge_reference_files:main',
        ],
    },
    include_package_data=True,
    package_data={
        '': ['resources/*', 'resources/omim/*', 'resources/candidate_lists/*', 'resources/references/*'],
    },
    license="MIT",
    description="A streamlit app for visualizing germline copy-number variants.",
    author="Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft",
    author_email="jerkrause@ukaachen.de",
    url="https://github.com/jerkrause/CNVizard",
)
