import numpy as np
import backend
from numpy import nan
import astropy.coordinates as c
import astropy.units as u

from format_unc_html import format_html

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, Range1d, HoverTool, Span, Band, Div, Spacer, Button, CustomJS 
from bokeh.plotting import figure
from bokeh.models import BoxZoomTool, ResetTool
import pickle
from format_unc_html import format_html

# Set up data
name = 'GRO J1744-28 @7..9kpc'
res0 = backend.get_bj_posterior(name)
dbins = backend.dbins
nh_conv_fac = backend.calib_nhw00_mean[0]

ebvl,ebv,ebvh = res0[1][0]
ebve = np.sqrt((ebvh-ebv)**2+(ebv-ebvl)**2)

avl,av,avh = res0[1][1]
ave = np.sqrt((avh-av)**2+(av-avl)**2)

akl,ak,akh = res0[1][2]
ake = np.sqrt((akh-ak)**2+(ak-akl)**2)

nhl,nhm,nhh = (res0[1][3].astype(np.float64))*1e21 # originally was in 1e21 units, we want 1e21, need to avoid overflows hence f64
nhe = 1e21*np.sqrt((nhh/1e21-nhm/1e21)**2+(nhm/1e21-nhl/1e21)**2)

dlow,dmean,dhi = res0[2:5]
gid = res0[6]

# h1 = res0[8]
# global source 
source = ColumnDataSource(data=dict(dist=dbins, 
                                    nh=nhm,
                                    nh21=nhm/1e21, 
                                    nh21e=nhe/1e21,
                                    nho21=ebv*nh_conv_fac, 
                                    nho21e=ebve*nh_conv_fac, 
                                    ebv=ebv,ebve=ebve, 
                                    av=av,ave=ave, 
                                    ak=ak,ake=ake,
                                    nhlb=np.maximum(np.ones_like(nhl)*1e20,nhl),nhub=nhh,
                                    nholb=np.maximum(np.ones_like(ebve)*1e20,ebvl*1e21*nh_conv_fac),nhoub=ebvh*1e21*nh_conv_fac))
                                    # h1_hi4pi=h1[0],h1_lab=h1[1],h1_dl=h1[2]))



distnce_span_med = Span(dimension='height',location=dmean, line_color='#B3A273', line_width=1.5, line_alpha=1.0)
distnce_span_lo = Span(dimension='height',location=dlow, line_color='#B3A273', line_width=1, line_alpha=0.5)
distnce_span_hi = Span(dimension='height',location=dhi, line_color='#B3A273', line_width=1, line_alpha=0.5)
# print(dlow,dhi)

plot = figure(height=550, width=600, title=name,
              tools=[BoxZoomTool(),ResetTool(),"save"],x_axis_label='Distance, pc', y_axis_label=r'\(N_H, \text{atoms/cm}^{-2}\)', x_axis_type="log", y_axis_type="log",x_range=(100,27000),y_range=(1e20,nhh.max()*2))

plot.line('dist', 'nh', source=source, line_width=3, line_alpha=0.5, line_color='#902535')

# plot.line('dist', 'h1_hi4pi', source=source, line_width=1, line_alpha=0.5, line_color='#B3A273')
# plot.line('dist', 'h1_lab', source=source, line_width=1, line_alpha=0.5, line_color='#B3A273')
# plot.line('dist', 'h1_dl', source=source, line_width=1, line_alpha=0.5, line_color='#B3A273')

band = Band(base="dist", lower="nhlb", upper="nhub", source=source,
            fill_alpha=0.1, fill_color="#B3A273", line_color="#B3A273")
band_optimistic = Band(base="dist", lower="nholb", upper="nhoub", source=source,
            fill_alpha=0.2, fill_color="#B3A273", line_color="#B3A273")

plot.add_layout(band)
plot.add_layout(band_optimistic)
plot.y_range.end = nhh.max()*2


plot.add_tools(HoverTool(
    tooltips = """
    <div style ="background-color:white;">         
            <div>
                <span style="font-size: 12px; color: #B3A273;">Distance: @dist{0.0} pc</span>
            </div>
            <div>
                <span style="font-size: 12px; color: #B3A273;"> E(B-V): @ebv±@ebve mag</span>
            </div>
            <div>
                <span style="font-size: 12px; color: #B3A273;"> A<sub>V</sub>: @av±@ave mag</span>
            </div>
            <div>
                <span style="font-size: 12px; color: #B3A273;"> A<sub>K</sub>: @ak±@ake mag</span>
            </div>
        
            <div>
                <span style="font-size: 12px; color: #B3A273;"> N<sub>H</sub>: @nh21±@nh21e &times 10<sup>21</sup>atoms cm<sup>-2</sup></span>
            </div>
            <div>
                <span style="font-size: 12px; color: #B3A273;"> N<sub>H,E(B-V)</sub>: @nho21±@nho21e &times 10<sup>21</sup>atoms cm<sup>-2</sup></span>
            </div>
        </div>
    """,
    
    mode='vline',
))

plot.add_layout(distnce_span_med)
plot.add_layout(distnce_span_lo)
plot.add_layout(distnce_span_hi)

#B3A273 - gold
#902535 - Karminrot

# Print header
header = Div(text="""
<h1>3DN<sub>H</sub>-tool: <span class="smaller-text">column density/extinction along line of sight for a given position</span>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<a href=https://uni-tuebingen.de><img src='http://astro.uni-tuebingen.de/~xrbcat/UT_WBMW_Rot_RGB_01.png' width=90px></a></h1>
<style>
  .smaller-text {
    font-size: 0.4em; /* Adjust the font size as needed */
    color: #902535; /* Optional: Change the color if desired */
  }
</style>
<p>This tool is a convinience addition to the pan-Galactic 3d exctinction/absorbtion maps presented in <a href=http://ads>Doroshenko 2024</a> aimed to enable quick look-up of extinction and equivalent absorbtion column for individual sources without downloading the <a href=http://zenodo>full map.</a> This work is based on and combines the results by <a href=https://ui.adsabs.harvard.edu/abs/2023arXiv230801295E/abstract>Edenhofer et al 2023</a> (up to 2 kpc), and a combination of <a href=https://ui.adsabs.harvard.edu/abs/2017ApJ...835...29Y/abstract>YMW16</a> and <a href=https://ui.adsabs.harvard.edu/abs/2016A%26A...596A.109P/abstract>Planck GNILC</a> maps to extend the maps in the Galactic plane beyond this distance, and then calibrated using several independent datasets. Please consider citing these publications if you use this tool. 
</p>
""", width=1000)
footer = Div(text = """
<div>
    <span style="font-size: 12px; color: #C9C9C9; text-align: right;"> 
by <a href=https://doroshv.github.io>Victor Doroshenko</a>, 2023, <a href=https://uni-tuebingen.de/index.php?id=3130>University of Tübingen</a>  <a href=https://uni-tuebingen.de/en/334>(Imprint)</a>, powered by <a href=https://bokeh.org> <img src=https://static.bokeh.org/logos/logotype.svg width=50></a>.</span></div>
""", width=1000, align='end')


def result_summary_text(res):
    """Generate summary for a source"""
    gc = res[7].transform_to(c.Galactic)
    log_of_nh = np.floor(np.log10(res[0][3][1]*1e21))
    tex = """<h3>Results:</h3><strong>Assumed position:&nbsp</strong><br> <em>(l,b):</em>%.2f&deg,%.2f&deg<br><em>(Ra,Dec): </em>%.2f&deg,%.2f&deg<br><em>Gaia ID:</em> %s<br>
                               %s 
                               %s %s %s %s %s"""%(gc.l.deg,gc.b.deg,res[7].ra.deg,res[7].dec.deg,res[6],
                              format_html(res[3],res[3]-res[2],res[4]-res[3],"<strong>Assumed distance:&nbsp</strong>"," pc<br>"),
                              format_html(res[0][3][1]*1e21/(10**log_of_nh),1e21*(res[0][3][1]-res[0][3][0])/(10**log_of_nh),1e21*(res[0][3][2]-res[0][3][1])/(10**log_of_nh),"<strong>Estimated N<sub>H</sub>:&nbsp </strong>","<div>&times 10<sup>%d</sup> cm<sup>-2</sup></div>"%log_of_nh),
                              format_html(res[0][0][1]*1e21*nh_conv_fac/(10**log_of_nh),1e21*nh_conv_fac*(res[0][0][1]-res[0][0][0])/(10**log_of_nh),1e21*nh_conv_fac*(res[0][0][2]-res[0][0][1])/(10**log_of_nh),"<strong>or N<sub>H,E(B-V)</sub>:&nbsp </strong>","<div>&times 10<sup>%d</sup> cm<sup>-2</sup></div>"%log_of_nh),
                              format_html(res[0][0][1],(res[0][0][1]-res[0][0][0]),(res[0][0][2]-res[0][0][1]),"<strong>Estimated E<sub>(B-V)</sub>:&nbsp </strong>"," mag"),
                              format_html(res[0][1][1],(res[0][1][1]-res[0][1][0]),(res[0][1][2]-res[0][1][1]),"<strong>Estimated A<sub>V</sub>:&nbsp </strong>"," mag"),
                              format_html(res[0][2][1],(res[0][2][1]-res[0][2][0]),(res[0][2][2]-res[0][2][1]),"<strong>Estimated A<sub>K</sub>:&nbsp </strong>"," mag"))
    return tex + """<br>
    <span class="smaller-text">The 1 uncertainties for N<sub>H</sub> E(B-V) and A<sub>V,K</sub> are conservative and account for estimated statistical and systematic uncertainties both in extinction and assumed distance range. 
    The N<sub>H,E(B-V)</sub> is calculated directly from E(B-V) and does not account for uncertainties in N<sub>H</sub>(E(B-V)) calibration.
    If distance is not given explicitly and <em>Gaia</em> counterpart within 2'' can be identified, approximated geometric priors by <a href=https://ui.adsabs.harvard.edu/abs/2021AJ....161..147B/abstract>Bailer-Jones et al 2021</a> are used. Check whether association with Gaia source is correct and you're fine with this assumption. Please refer to the curve on the right if in doubt (vertical lines indicate distance used and band estimated absorbtion column uncertainties).</span>
    <style>
      .smaller-text {
        font-size: 0.7em; /* Adjust the font size as needed */
        color: #C9C9C9; /* Optional: Change the color if desired */
      }
    </style>
    
    """
    
result_summary = Div(text=result_summary_text(res0), width=300)

# styles=[".b_export button { background-color: #B3A273 !important; }"]
download_button = Button(label="Export curve (CSV)", button_type="default",width=100)
download_button.js_on_event("button_click", CustomJS(args=dict(source=source),
                            code=open("nh3d/download.js").read()))
                            

spacer = Spacer(width=100)

# Set up widgets
text_explanation = Div(text="""<h3>Enter sky coordinates or source name for the position of interest:<br>
 <span class="smaller-text">needs to be resolvable by <em>Simbad</em>, i.e. "Cen X-3", "170.31 -60.62", "gal 0 -0.5 @3..5kpc" etc.</span></h3>
<style>
  .smaller-text {
    font-size: 0.7em; /* Adjust the font size as needed */
    color: #C9C9C9; /* Optional: Change the color if desired */
  }
</style>

""",width=550)
text = TextInput(title="", value='GRO J1744-28 @8kpc', width=300)

# Set up callbacks
def update_title(attrname, old, new):
    global source
    source.data = dict(dist=dbins, nh=np.nan*dbins,nh21=np.nan*dbins, nh21e=np.nan*dbins, ebv=np.nan*dbins,ebve=np.nan*dbins, av=np.nan*dbins,ave=np.nan*dbins, ak=np.nan*dbins,ake=np.nan*dbins,nhlb=np.nan*dbins,nhub=np.nan*dbins)
    result_summary.text = "<h3>Can not resolve position, check input</h3>"
    plot.visible=False
    res0 = backend.get_bj_posterior(text.value)
    if res0:
        plot.title.text = text.value
        ebvl,ebv,ebvh = res0[1][0]
        ebve = np.sqrt((ebvh-ebv)**2+(ebv-ebvl)**2)

        avl,av,avh = res0[1][1]
        ave = np.sqrt((avh-av)**2+(av-avl)**2)

        akl,ak,akh = res0[1][2]
        ake = np.sqrt((akh-ak)**2+(ak-akl)**2)

        nhl,nhm,nhh = (res0[1][3].astype(np.float32))*1e21 # originally was in 1e21 *units*, we want true 1e21
        nh = nhm
        nhe = 1e21*np.sqrt((nhh/1e21-nhm/1e21)**2+(nhm/1e21-nhl/1e21)**2)

        dlow,dmean,dhi = res0[2:5]
        gid = res0[6]
        # h1 = res0[8]
        source.data=dict(dist=dbins, 
                                    nh=nhm,
                                    nh21=nhm/1e21, 
                                    nh21e=nhe/1e21,
                                    nho21=ebv*nh_conv_fac, 
                                    nho21e=ebve*nh_conv_fac, 
                                    ebv=ebv,ebve=ebve, 
                                    av=av,ave=ave, 
                                    ak=ak,ake=ake,
                                    nhlb=np.maximum(np.ones_like(nhl)*1e20,nhl),nhub=nhh,
                                    nholb=np.maximum(np.ones_like(ebve)*1e20,ebvl*1e21*nh_conv_fac),nhoub=ebvh*1e21*nh_conv_fac)
                                    # h1_hi4pi=h1[0],h1_lab=h1[1],h1_dl=h1[2])
        distnce_span_med.location=dmean
        distnce_span_lo.location=dlow
        distnce_span_hi.location=dhi
        plot.y_range.start = 1e20
        plot.y_range.end = nhh.max()*2
        result_summary.text = result_summary_text(res0)
        plot.visible=True
        # plot.tools[1].update()
        download_button.visible=True
    else:
        plot.visible=False
        download_button.visible=False
        result_summary.text = "<h3>Can not resolve position, please check your input</h3> allowed formats are Simbad-resolvable identifiers, equatorial coordinates in degrees and hms/dms form (i.e. 16 17 04.4 -15 31 14.83 or 16h17m04.4s -15d31m14.83s), Galactic coordinates in degrees prepended by `gal', and everything that astropy can parse. To indicate distance range for which you would like to estimate absorbtion column use syntax @5..8kpc or @100pc, i.e. indicate units and use .. to specify ranges."
        plot.title.text = "Can not resolve position"
        y =  np.array([np.nan])
        ye = np.array([np.nan])
        nh = y
        nhe = ye
        dlow,dmean,dhi = np.array([np.nan]*3)
        gid = 'not found'
        nhlb = nh-nhe
        source.data = dict(dist=dbins, nh=nh,nh21=nh/1e21, nh21e=nhe/1e21, av=y, ave=ye,nhlb=nhlb,nhub=nh+nhe)
        

text.on_change('value', update_title)

# Set up layouts and add to document
inputs = row(text_explanation,column(Spacer(height=15),text))
curdoc().add_root(row(spacer,column(Spacer(height=20),header,Spacer(height=10),column(inputs, Spacer(width=10), row(column(result_summary,download_button),Spacer(width=10),column(Spacer(height=10),plot)),Spacer(height=20),footer), width=900)))
