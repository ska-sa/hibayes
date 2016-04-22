import sys
import pprint

try:
    sys.path.append("/workdata/leda/hibayes")
except:
    pass

from hibayes import parse_config


def test_parse_config():
    cd = parse_config.parse_config("test_config.ini")
    pprint.pprint(cd)

if __name__ == "__main__":
    test_parse_config()