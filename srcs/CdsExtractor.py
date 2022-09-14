import Bio.Seq
import Bio.SeqFeature


def extract_cds(seq: Bio.Seq.Seq, feature_location: Bio.SeqFeature.FeatureLocation):
    return feature_location.extract(seq)


if __name__ == '__main__':
    pass
