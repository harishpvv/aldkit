from setuptools import setup, find_packages
import glob
  
with open("README.md", 'r') as f:
    long_description = f.read()
  
setup(
        name ='aldkit',
        version ='1.1.10',
        author ='Harish PVV',
        author_email ='harishpvv@gmail.com',
        description ="Parser for anharmonic lattice dynamics package developed by Prof. Ankit Jain.",
        long_description = long_description,
        long_description_content_type ="text/markdown",
        license ='MIT',
        packages = find_packages(),
        entry_points ={
            'console_scripts': [
                'ald_input = aldkit.ald_input:main',
                'ald_pbs = aldkit.ald_pbs:main',
                'ald_help = aldkit.ald_help:main'
            ]
        },
        classifiers =(
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ),
        keywords ='ald',
        install_requires = ['ase', 'spglib'],
        zip_safe = False,
        data_files = [("pseudos/PBE_ONCV/", glob.glob("src/aldkit/pseudos/PBE_ONCV/*.UPF")),
                      ("pseudos/LDA_ONCV/", glob.glob("src/aldkit/pseudos/LDA_ONCV/*.UPF")),
                      ("pseudos/PBESOL_ONCV/", glob.glob("src/aldkit/pseudos/PBESOL_ONCV/*.UPF")),
                      ("aldcodes/", glob.glob("src/aldkit/aldcodes/*"))],
        include_package_data = True

        )
