import pickle,sys,warnings
import numpy as np
import matplotlib.pylab as plt
from math import pi,degrees
import PySimpleGUI as sg
from numpy import sin,cos
from numpy.linalg import eig, inv
from matplotlib.patches import Circle,Ellipse
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
plt.rcParams.update({'text.usetex':False})
warnings.filterwarnings("ignore")
#Ellipse functions from https://stackoverflow.com/questions/13635528/fit-a-ellipse-in-python-given-a-set-of-points-xi-xi-yi
def EllipseFit(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    return( V[:,n] )
def EllipseCenter(ellipse):
    b,c,d,f,g,a = ellipse[1]/2,ellipse[2],ellipse[3]/2,ellipse[4]/2,ellipse[5],ellipse[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])
def EllipseAngle(ellipse):
    b,c,d,f,g,a = ellipse[1]/2,ellipse[2],ellipse[3]/2,ellipse[4]/2,ellipse[5],ellipse[0]
    return 0.5*np.arctan(2*b/(a-c))
def EllipseAxes(ellipse):
    b,c,d,f,g,a = ellipse[1]/2,ellipse[2],ellipse[3]/2,ellipse[4]/2,ellipse[5],ellipse[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1,res2])
#Masking Functions
def MaskCircle(rad,cen_x,cen_y,isophote,mode):
    assert mode in ['Inclusive','Exclusive'], 'Masking Mode Error'
    masked_iso = [[],[]]
    for i in np.arange(len(isophote[0])):
        r = np.sqrt((isophote[0][i]-500-cen_y)**2+(isophote[1][i]-500-cen_x)**2)
        if r<rad and mode=='Inclusive':
            masked_iso[0].append(isophote[0][i])
            masked_iso[1].append(isophote[1][i])
        if r>rad and mode=='Exclusive':
            masked_iso[0].append(isophote[0][i])
            masked_iso[1].append(isophote[1][i])
    masked_iso = ( np.array(masked_iso[0]),np.array(masked_iso[1]) )
    return(masked_iso)
def MaskSlice(minang,maxang,isophote,mode):
    ###NOT READY, angle shit
    assert mode in ['Inclusive','Exclusive'], 'Masking Mode Error'
    masked_iso = [[],[]]
    for i in np.arange(len(isophote[0])):
        y,x = isophote[0][i]-500,isophote[1][i]-500
        if x<0:
            if y>0: ang = degrees(np.arctan(y/x))+180
            if y<0: ang = degrees(np.arctan(y/x))-180
        else: ang = degrees(np.arctan(y/x))
        if minang<ang<maxang and mode=='Inclusive':
            masked_iso[0].append(isophote[0][i])
            masked_iso[1].append(isophote[1][i])
        if (minang>ang or maxang<ang) and mode=='Exclusive':
            masked_iso[0].append(isophote[0][i])
            masked_iso[1].append(isophote[1][i])
    masked_iso = ( np.array(masked_iso[0]),np.array(masked_iso[1]) )
    return masked_iso
def InitializeGUI(PlotName):
    #GUI Properties
    _VARS = {'window':False,'fig_agg':False,'pltFig':False}#,'a':np.nan,'b':np.nan}
    plt.style.use('Solarize_Light2')
    GuiFont = 'Any 16'
    GuiColor = '#E8E8E8'
    sg.theme('black')
    layout = [
        [sg.Text(text=PlotName,
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black')],
        [sg.Canvas(key='figCanvas',
            background_color=GuiColor)],
        [sg.Text(text='Circular Mask:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.Listbox(['Inclusive','Exclusive'],
            key='CMode'),
        sg.Text(text='Radius:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.Slider(range=(0,400), 
            orientation='h',size=(10,10),
            default_value=0,
            background_color=GuiColor,
            resolution=1,
            text_color='Black',
            key='RadiusAdjust',
            enable_events=True),
        sg.Text(text='Cen-X:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.Slider(range=(-150,150), 
            orientation='h',size=(10,10),
            default_value=0,
            background_color=GuiColor,
            resolution=1,
            text_color='Black',
            key='CenXAdjust',
            enable_events=True),
        sg.Text(text='Cen-Y:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.Slider(range=(-150,150), 
            orientation='h',size=(10,10),
            default_value=0,
            background_color=GuiColor,
            resolution=1,
            text_color='Black',
            key='CenYAdjust',
            enable_events=True)],
        [sg.Text(text='Angular Mask:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.Listbox(['Exclusive','Inclusive'],
            key='AMode'),
        sg.Text(text='Min Angle:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.Slider(range=(-180,180), 
            orientation='h',size=(17,10),
            default_value=0,
            background_color=GuiColor,
            resolution=1,
            text_color='Black',
            key='MinAngle',
            enable_events=True),
        sg.Text(text='Max Angle:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.Slider(range=(-180,180), 
            orientation='h',size=(17,10),
            default_value=0,
            background_color=GuiColor,
            resolution=1,
            text_color='Black',
            key='MaxAngle',
            enable_events=True)],
        [sg.Button('Exit',font=GuiFont),
        sg.Button('Ignore',font=GuiFont,pad=((0,130),(0,0))),
        sg.Text(text='Isophote %:',
            font=GuiFont,
            background_color=GuiColor,
            text_color='Black'),
        sg.InputText('1',
            background_color=GuiColor,
            font=GuiFont,
            size=3,
            text_color='Black',
            key='Iso%'),
        sg.Button('Reset',font=GuiFont),
        sg.Button('Mask', font=GuiFont),
        sg.Button('Fit Ellipse', font=GuiFont),
        sg.Button('Save', font=GuiFont)]
    ]
    _VARS['window'] = sg.Window('Isophote Fitting',
                                    layout,
                                    finalize=True,
                                    resizable=True,
                                    location=(100, 100),
                                    element_justification="center",
                                    background_color=GuiColor)
    return _VARS
def draw_figure(canvas, figure):
        figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
        return figure_canvas_agg
def drawChart():
    _VARS['pltFig'] = plt.figure()
    plt.grid(b=None)
    _VARS['pltFig'].axes[0].set_xlim([0,1000])
    _VARS['pltFig'].axes[0].set_ylim([0,1000])
    plt.imshow(im)
    plt.scatter(500,500,marker='+',s=10**2,c='w')
    plt.scatter(iso[1],iso[0],c='r',s=.5**2)
    _VARS['fig_agg'] = draw_figure(_VARS['window']['figCanvas'].TKCanvas,_VARS['pltFig'])
def updateChart(DrawCircle=False,DrawAngle=False,DrawEllipse=False):
    _VARS['fig_agg'].get_tk_widget().forget()
    plt.clf()
    plt.grid(b=None)
    _VARS['pltFig'].axes[0].set_xlim([0,1000])
    _VARS['pltFig'].axes[0].set_ylim([0,1000])
    plt.imshow(im)
    plt.scatter(500,500,marker='+',s=10**2,c='w')
    plt.scatter(iso[1],iso[0],c='r',s=.5**2)
    if DrawCircle:
        _VARS['pltFig'].axes[0].add_patch(Circle((500+values['CenXAdjust'],500+values['CenYAdjust']),
                                        values['RadiusAdjust'],color='w',fill=False))
    if DrawAngle:
        plt.plot([500,500+710*np.cos(np.radians(values['MinAngle']))],
                 [500,500+710*np.sin(np.radians(values['MinAngle']))],c='w',linewidth=1)
        plt.plot([500,500+710*np.cos(np.radians(values['MaxAngle']))],
                 [500,500+710*np.sin(np.radians(values['MaxAngle']))],c='w',linewidth=1)
    if DrawEllipse:
        E = EllipseFit(iso[1],iso[0])
        cen = EllipseCenter(E)
        phi = EllipseAngle(E)
        a,b = EllipseAxes(E)
        #Plot the ellipse fit on the image and set image title to axis ratio
        _VARS['pltFig'].axes[0].add_patch(Ellipse(cen,2*a,2*b,angle=degrees(phi),facecolor='None',edgecolor='orange'))
        plt.plot([-a*cos(phi)+cen[0],a*cos(phi)+cen[0]],[-a*sin(phi)+cen[1],a*sin(phi)+cen[1]],
                color='orange')
        plt.plot([-b*cos(phi+pi/2)+cen[0],b*cos(phi+pi/2)+cen[0]],[-b*sin(phi+pi/2)+cen[1],
                b*sin(phi+pi/2)+cen[1]],color='orange')   
        Data[halo][view] =  min([a,b])/max([a,b])
    _VARS['fig_agg'] = draw_figure(_VARS['window']['figCanvas'].TKCanvas,_VARS['pltFig'])


Images = pickle.load(open('../DataFiles/ObservedShapes.ImageData.pickle','rb'))
outfile = '../DataFiles/ObservedShapes.ShapeData.pickle'
Data = {}

for halo in Images:
    Data[halo] = {}
    for view in ['Faceon','Sideon']:
        Data[halo][view] = np.nan
        v,im = Images[halo][f'SB_{view}'],Images[halo][f'Im_{view}']
        tol = .01
        iso = np.where((im>v*(1-tol)) & (im<v*(1+tol)))

        DC,DA,DE = False,False,False
        _VARS = InitializeGUI(f'Halo {halo} - {view}')
        drawChart()
        while True:
            event, values = _VARS['window'].read()
            if event in [sg.WIN_CLOSED,'Exit']:
                print('Aborting Code')
                sys.exit(0)
            if event in ['Ignore','Save']:
                break
            if event=='Reset':
                tol = float(values['Iso%'])/100
                iso = np.where((im>v*(1-tol)) & (im<v*(1+tol)))
                DC,DA,DE = False,False,False
                updateChart(DrawCircle=DC,DrawAngle=DA,DrawEllipse=DE)
            if event=='Mask':
                DE = False
                if DC:
                    iso = MaskCircle(values['RadiusAdjust'],cen_x=values['CenXAdjust'],
                            cen_y=values['CenYAdjust'],isophote=iso,mode=values['CMode'][0])
                if DA:
                    iso = MaskSlice(values['MinAngle'],values['MaxAngle'],
                            isophote=iso,mode=values['AMode'][0])
                updateChart(DrawCircle=DC,DrawAngle=DA,DrawEllipse=DE)
            if event in ['RadiusAdjust','CenXAdjust','CenYAdjust'] and values['CMode']!=[]:
                DC,DE = True,False
                updateChart(DrawCircle=DC,DrawAngle=DA,DrawEllipse=DE)
            if event in ['MinAngle','MaxAngle'] and values['AMode']!=[]:
                DA,DE = True,False
                updateChart(DrawCircle=DC,DrawAngle=DA,DrawEllipse=DE)
            if event=='Fit Ellipse':
                DE = True
                updateChart(DrawCircle=DC,DrawAngle=DA,DrawEllipse=True)
        plt.close()
        _VARS['window'].close()

#out = open(outfile,'wb')
#pickle.dump(Data,out)
#out.close()
#print('File Saved.')
print(Data)