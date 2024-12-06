from hapi.utils.reader_of_yaml import YamlReader


def test_yamlreader():
    reader = YamlReader("./testfiles/testconfig.yaml")
    position_list_reference = reader.position_list_reference
    position_list_deletion = reader.position_list_deletion
    top4_snps_list = reader.top4_snps_list
    chrom = reader.chromosome

    assert position_list_reference == [
        [46414944, 46414975],
        [46414945, 46414976],
        [46414946, 46414977],
        [46414947, 46414978],
    ]
    assert position_list_deletion == [46414943]
    assert top4_snps_list == ["rs113341849", "rs113010081", "rs11574435", "rs79815064"]
    assert chrom == str(3)
