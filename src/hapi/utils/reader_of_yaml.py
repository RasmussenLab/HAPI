from typing import List


class YamlReader:
    def __init__(self, file_path: str) -> None:
        import yaml

        with open(file_path, "r") as f:
            self.yaml_loader = yaml.safe_load(f)

    @property
    def position_list_reference(self) -> List[List[int]]:
        return self.yaml_loader["position_list_reference"]

    @property
    def position_list_deletion(self) -> list:
        return [self.yaml_loader["position_list_deletion"]]

    @property
    def top4_snps_list(self) -> list:
        return self.yaml_loader["top4_snps_list"]

    @property
    def chromosome(self) -> str:
        return str(self.yaml_loader["chromosome"])

    @property
    def deletion_length(self) -> int:
        return int(self.yaml_loader["deletion_length"])
