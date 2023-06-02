import os
from pathlib import Path


class AssetsManager:
    def __init__(self, assets_dir=None):
        self.assets_dir = os.path.dirname(os.path.abspath(__file__)) if assets_dir is None else assets_dir
        self._load_assets()

    @property
    def all_assets(self):
        return [
            attr for attr in dir(self)
            if not attr.startswith('_')
            and attr not in ['assets_dir', '_load_assets', 'all_assets']
        ]

    def _load_assets(self):
        # TODO: Replace with os.scandir() for Python 3.5+
        for root, dirs, files in os.walk(self.assets_dir):
            for file in files:
                if file.endswith('.txt') and not file.startswith('._'):
                    attr_name = file[:-4]  # remove .txt extension
                    attr_name = attr_name.replace('-', '_')  # replace hyphen with underscore for valid attribute name
                    attr_path = Path(root) / file
                    with open(attr_path, 'r') as f:
                        try:
                            setattr(self, attr_name, [line.strip() for line in f if line.strip()])
                        except UnicodeDecodeError:
                            print(f'UnicodeDecodeError: {attr_path}')

data_assets = AssetsManager()