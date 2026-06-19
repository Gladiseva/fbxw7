import argparse
from pathlib import Path

try:
    import yaml
except ImportError as exc:
    raise SystemExit(
        "PyYAML is required to read the config file. "
        "Install the manuscript environment with `conda env create -f environment.yml`."
    ) from exc


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = REPO_ROOT / "config" / "default.yaml"


def parse_config_arg(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG_PATH),
        help="Path to a YAML configuration file.",
    )
    return parser.parse_args()


def load_config(config_path):
    config_path = Path(config_path).expanduser()
    if not config_path.is_absolute():
        config_path = (Path.cwd() / config_path).resolve()

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    project_root = Path(config.get("project_root", ".")).expanduser()
    if not project_root.is_absolute():
        project_root = (config_path.parent / project_root).resolve()

    config["_project_root"] = project_root
    config["_config_path"] = config_path
    return config


def resolve_path(config, key):
    value = config["paths"][key]
    path = Path(value).expanduser()
    if path.is_absolute():
        return str(path)
    return str((config["_project_root"] / path).resolve())


def get_markers(config):
    return list(config["analysis"]["markers"])
