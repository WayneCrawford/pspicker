import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()
    
version={}
with open("ps_picker/version.py") as fp:
    exec(fp.read(),version)

setuptools.setup(
    name="ps_picker",
    version=version['__version__'],
    author="Christian Baillard",
    author_email="crawford@ipgp.fr",
    description="P and S wave picker using kurtosis",
    long_description=long_description,
    long_description_content_type="text/x-rst; charset=UTF-8",
    url="https://github.com/WayneCrawford/ps_picker",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
          'obspy>=1.2',
          'numpy',
          'scipy',
          'matplotlib',
          'pyyaml>=3.0,<4',
          'jsonschema>=2.6',
          'jsonref>=0.2'
      ],
    entry_points={},
    python_requires='>=3.6',
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    classifiers=(
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ),
    keywords='seismology'
)
