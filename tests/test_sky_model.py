import sys
import pprint

try:
    sys.path.append("/workdata/leda/hibayes")
except:
    pass

from hibayes import parse_config
from hibayes import sky_model

def test_sky_model():
    rp = parse_config.parse_config("test_config.ini")

    Tmeas, freqs = sky_model.generate_simulated_data(rp, plot_data=True)

if __name__ == "__main__":
    test_sky_model()