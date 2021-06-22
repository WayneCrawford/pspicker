import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    
version={}
with open("pspicker/version.py") as fp:
    exec(fp.read(),version)

setuptools.setup(
    name="pspicker",
    version=version['__version__'],
    author="Christian Baillard",
    author_email="crawford@ipgp.fr",
    description="Kurtosis-based P and S wave picker",
    long_description=long_description,
    long_description_content_type="text/markdown; charset=UTF-8",
    url="https://github.com/WayneCrawford/pspicker",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
          'obspy>=1.2',
          'numpy>=1.18',
          'scipy>=1.5',
          'verboselogs>=1.7',
          'matplotlib>=3.2',
          'pyyaml>=3.0',
          'jsonschema>=2.6',
          'jsonref>=0.2'
      ],
    entry_points={},
    python_requires='>=3.8',
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
