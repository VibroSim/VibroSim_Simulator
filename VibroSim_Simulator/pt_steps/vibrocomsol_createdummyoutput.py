import os
import os.path

try:
    # py2.x
    from urllib import pathname2url
    from urllib import url2pathname
    from urllib import quote
    from urllib import unquote
    pass
except ImportError:
    # py3.x
    from urllib.request import pathname2url
    from urllib.request import url2pathname
    from urllib.parse import quote
    from urllib.parse import unquote
    pass

from limatix.dc_value import hrefvalue as hrefv

def run(dc_dest_href,dc_measident_str):

    dummy_heatingdata_href = hrefv(quote(dc_measident_str+"_dummy_heatingdata.txt"),dc_dest_href)
    
    dummy_heatingdatapath = dummy_heatingdata_href.getpath()

    # Create directory if needed
    if not os.path.exists(os.path.split(dummy_heatingdatapath)[0]):
        os.mkdir(os.path.split(dummy_heatingdatapath)[0])
        pass

    fh = open(dummy_heatingdatapath,"w")
    fh.write("% t(s) \t r(m) \t side1_heating(W/m^2) \t side2_heating(W/m^2)\n")
    fh.write("0\t0\t0\t0\n")
    fh.close()
    return { "dc:dummy_heatingdata": dummy_heatingdata_href }
