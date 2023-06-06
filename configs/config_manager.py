import json
import os


# TODO: Deprecate in JSON reading of plotting config paths

class ConfigManager:
    def __init__(self, base_path=None, config_path='config.json'):
        self.base_path = os.path.dirname(os.path.abspath(__file__)) if base_path is None else base_path

        # Load config.json
        with open(os.path.join(self.base_path, config_path)) as f:
            self.config = json.load(f)

        # Convert paths to absolute paths
        for key in ['data', 'plotting', 'plots']:
            setattr(self, key, self._load_data_configs(key))

    def _load_data_configs(self, key):
        return {k: self._to_absolute_path(v) for k, v in self.config[key].items()}

    def _to_absolute_path(self, config_item):
        if isinstance(config_item, dict):
            return {k: self._to_absolute_path(v) for k, v in config_item.items()}

        elif isinstance(config_item, str) and config_item.startswith('./'):
            proj_base_path = os.path.dirname(self.base_path)
            return os.path.join(proj_base_path, config_item[2:])

        else:
            return config_item

config = ConfigManager()