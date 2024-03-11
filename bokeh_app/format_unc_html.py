# python -m pip install uncertainties
from uncertainties import ufloat
import re
from numpy import *

value, low_err, upper_err = [1,2.5,3.02,4.005,4.1],[1.1,0.3,0.02,0.05,0.03],[0.3,1.2,0.05,0.01,0.03]

def format_html(value,low_err,upper_err,prefix,suffix):
    """Format value and uncertainties to latex form. It makes sense to adjust for exponent before formatting (i.e. divide by 1e-12 or whatever ant then modify output accordingly)"""
    css = lambda val,le,ue,prefix,suffix: """<style>.expression{
  display:flex;
  align-items: center;
}
.supsub {
  display: flex;
  flex-direction: column;
  margin-left: 2px;
  margin-right: 2px;
  }
.subscript {
  color: black; 
  display:flex;
  }
  
.superscript {
  color: black; 
  display:flex; 
  }
    </style>
  <div class="expression">
  %s %s 
  <span class='supsub'>
    <sup class='superscript'>%s</sup>
    <sub class='subscript'>%s</sub>
  </span> %s
  </div>

  """%(prefix,val,ue,le,suffix)
    value, low_err, upper_err = array(value,dtype=float,ndmin=1), array(low_err,dtype=float,ndmin=1),array(upper_err,dtype=float,ndmin=1)
    if low_err+upper_err==0:
        return "<div> %s %.1f %s</div>"%(prefix,value,suffix)
    str_lowerr = [ufloat(x).format("%f").replace('+/-','<sub>-')+'</sub>' for x in zip(value,low_err)]
    str_uperr = [ufloat(x).format("%f").replace('+/-','<sup>+')+'</sup>' for x in zip(value,upper_err)]
    formatted = [x[0]+re.findall('<sup>.*</sup>',x[1])[0] for x in zip(str_lowerr,str_uperr)]
    for i in range(len(formatted)):
        k = formatted[i]
        vv =  k[:k.find('<')]
        perr = re.findall('<sup>\+(.*)</sup>',k)[0]
        merr = re.findall('<sub>\-(.*)</sub>',k)[0]
        if merr==perr:
                formatted[i]=vv+r"(%s)"%merr
    formatted = formatted[0]
    try:
        v = formatted[:formatted.find('<')]
        le = re.findall('<sub>(.*)</sub>',formatted)[0]
        ue = re.findall('<sup>(.*)</sup>',formatted)[0]
        formatted = css(v,le,ue,prefix,suffix)
    except:
        formatted = "<div>%s %s %s</div>"%(prefix,formatted,suffix)
    return formatted

