import venv
import os
import sys
import subprocess


class TLSEnvBuilder(venv.EnvBuilder):
    """
    This builder creates a virtual Python environment in a specified directory,
    and also sets up R with specified packages.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def post_setup(self, context):
        os.environ['VIRTUAL_ENV'] = context.env_dir
        self.setup_r(context.env_dir)
        print("Virtual environment and R setup completed successfully.")

    def setup_r(self, env_dir):
        # Path to the requirements.R file
        requirements_r = os.path.join(os.path.dirname(__file__), 'requirements.R')

        # Create an R script to install packages
        script_path = os.path.join(env_dir, 'install_r_packages.R')
        with open(script_path, 'w') as f:
            f.write('options(repos = c(CRAN = "https://cloud.r-project.org"))\n')  # select a CRAN mirror
            with open(requirements_r) as reqs:
                for line in reqs:
                    package, source = line.strip().split(',')
                    if source == 'CRAN':
                        f.write(f'install.packages("{package}")\n')
                    elif source == 'Bioconductor':
                        f.write(f'if (!require("BiocManager", quietly = TRUE)) '
                                'install.packages("BiocManager")\n'
                                f'BiocManager::install("{package}")\n')

        # Execute the R script
        result = subprocess.run(["Rscript", script_path], check=True)
        if result.returncode != 0:
            raise Exception(f"Error in R setup: {result.stderr}")


def main(args=None):
    assert sys.version_info >= (3, 3), "This script is only for use with Python 3.3 or later"

    import argparse

    parser = argparse.ArgumentParser(prog=__name__,
                                     description='Creates virtual Python environments in one or '
                                                 'more target directories.')
    parser.add_argument('dirs', metavar='ENV_DIR', nargs='+',
                        help='A directory to create the environment in.')
    parser.add_argument('--clear', default=False, action='store_true',
                        dest='clear', help='Delete the contents of the '
                                           'environment directory if it '
                                           'already exists, before '
                                           'environment creation.')
    options = parser.parse_args(args)
    builder = TLSEnvBuilder(clear=options.clear)
    for d in options.dirs:
        builder.create(d)


if __name__ == '__main__':
    try:
        main()
        sys.exit(0)
    except Exception as e:
        sys.exit(f'Error: {e}')