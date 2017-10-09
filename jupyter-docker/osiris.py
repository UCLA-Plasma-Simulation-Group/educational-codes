import os
import shutil
import subprocess
import IPython.display
import h5py
import numpy as np
import matplotlib.pyplot as plt
#from os import listdir
from ipywidgets import interact
# %matplotlib inline

def runosiris(rundir='',inputfile='osiris-input.txt'):

    def execute(cmd):
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    workdir = os.getcwd()
    workdir += '/' + rundir
    print(workdir)
    if(not os.path.isdir(workdir)):
       os.mkdir(workdir)
    if(rundir != ''):
        shutil.copyfile('osiris-1D.e',workdir+'/osiris-1D.e')
        shutil.copyfile(inputfile,workdir+'/osiris-input.txt')
    for path in execute(["./osiris-1D.e","-w",workdir,"osiris-input.txt"]):
        IPython.display.clear_output(wait=True)
        print(path, end='')
    return


def field(rundir='',dataset='e1',time=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[]):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'FLD', dataset)
    files = os.listdir(odir)

    i = 0
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        if(fhere.attrs['TIME'] >= time):
            i = j
            break

    fhere = h5py.File(os.path.join(odir,files[i]), 'r')

    plt.figure(figsize=(6, 3.2))
    plt.title(dataset+' at t = '+str(fhere.attrs['TIME']))
    plt.xlabel('$x_1 [c/\omega_p]$')
    plt.ylabel(dataset)

    xaxismin = fhere['AXIS']['AXIS1'][0]
    xaxismax = fhere['AXIS']['AXIS1'][1]

    nx = len(fhere[dataset][:])
    dx = (xaxismax-xaxismin)/nx

    plt.plot(np.arange(0,xaxismax,dx),fhere[dataset][:])

    if(xlim != [-1,-1]):
        plt.xlim(xlim)
    if(ylim != [-1,-1]):
        plt.ylim(ylim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)

    plt.show()


def phasespace(rundir='',dataset='p1x1',species='electrons',time=0,
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[]):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'PHA', dataset, species)
    files = os.listdir(odir)

    i = 0
    for j in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[j]), 'r')
        if(fhere.attrs['TIME'] >= time):
            i = j
            break

    fhere = h5py.File(os.path.join(odir,files[i]), 'r')

    plt.figure(figsize=(6, 3.2))
    plt.title(dataset+' phasespace at t = '+str(fhere.attrs['TIME']))
    plt.xlabel('$x_1 [c/\omega_p]$')
    plt.ylabel('$p_1 [m_ec]$')

    if(len(fhere['AXIS']) == 1):

        xaxismin = fhere['AXIS']['AXIS1'][0]
        xaxismax = fhere['AXIS']['AXIS1'][1]

        nx = len(fhere[dataset][:])
        dx = (xaxismax-xaxismin)/nx

        plt.plot(np.arange(0,xaxismax,dx),fhere[dataset][:])

    elif(len(fhere['AXIS']) == 2):

        xaxismin = fhere['AXIS']['AXIS1'][0]
        xaxismax = fhere['AXIS']['AXIS1'][1]
        yaxismin = fhere['AXIS']['AXIS2'][0]
        yaxismax = fhere['AXIS']['AXIS2'][1]

        plt.imshow(np.log(np.abs(fhere[dataset][:,:]+1e-12)),
                   aspect='auto',
                   extent=[xaxismin, xaxismax, yaxismin, yaxismax])
        plt.colorbar(orientation='vertical')


    if(xlim != [-1,-1]):
        plt.xlim(xlim)
    if(ylim != [-1,-1]):
        plt.ylim(ylim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)

    plt.show()


def lineinteract(rundir='',dataset='e1',
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[]):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'FLD', dataset)
    files = os.listdir(odir)

    f0 = h5py.File(os.path.join(odir,files[0]), 'r')
    xaxismin = f0['AXIS']['AXIS1'][0]
    xaxismax = f0['AXIS']['AXIS1'][1]
    nx = len(f0[dataset][:])
    dx = (xaxismax-xaxismin)/nx
    xaxis = np.arange(0,xaxismax,dx)

    data = []
    for i in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[i]), 'r')
        data.append([fhere[dataset][:],fhere.attrs['TIME'],fhere.attrs['DT']])

    def fu(n):
        plt.figure(figsize=(8, 4))
        plt.plot(xaxis,data[n][0])
        plt.title('time = '+str(data[n][1]))
        if(xlim != [-1,-1]):
            plt.xlim(xlim)
        if(ylim != [-1,-1]):
            plt.ylim(ylim)
        if(zlim != [-1,-1]):
            plt.clim(zlim)
        return plt

    interact(fu,n=(0,len(data)-1))


def phaseinteract(rundir='',dataset='p1x1',species='electrons',
    xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[]):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'PHA', dataset, species)
    files = os.listdir(odir)

    data = []
    for i in range(len(files)):
        fhere = h5py.File(os.path.join(odir,files[i]), 'r')
        data.append([fhere[dataset][:,:],fhere.attrs['TIME'],fhere.attrs['DT']])
        xaxis = fhere['AXIS/AXIS1'][:]
        yaxis = fhere['AXIS/AXIS2'][:]

    def fu(n):
        plt.figure(figsize=(8, 4))
        plt.imshow(data[n][0],
               extent=[xaxis[0], xaxis[1], yaxis[0], yaxis[1]],
               aspect='auto')
        plt.title('time = '+str(data[n][1]))
        plt.colorbar()
        if(xlim != [-1,-1]):
            plt.xlim(xlim)
        if(ylim != [-1,-1]):
            plt.ylim(ylim)
        if(zlim != [-1,-1]):
            plt.clim(zlim)
        return plt

    interact(fu,n=(0,len(data)-1))


def xtplot(rundir='',dataset='e3',xlim=[-1,-1],ylim=[-1,-1],zlim=[-1,-1],
    plotdata=[]):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'FLD', dataset)
    files = os.listdir(odir)

    fhere = h5py.File(os.path.join(odir, files[0]), 'r')
    xaxismin = fhere['AXIS']['AXIS1'][0]
    xaxismax = fhere['AXIS']['AXIS1'][1]
    taxismin = fhere.attrs['TIME']

    fhere = h5py.File(os.path.join(odir, files[len(files)-1]), 'r')
    taxismax = fhere.attrs['TIME']

    nx = len(fhere[dataset][:])
    dx = (xaxismax-xaxismin)/nx
    nt = len(files)

    # fhere = h5py.File(odir+files[len(files)-1], 'r')
    fhere2 = h5py.File(os.path.join(odir, files[len(files)-2]), 'r')
    dt = fhere.attrs['TIME']-fhere2.attrs['TIME']
    # print(dt)

    if(plotdata == []):
        data = []
        for f in files:
            fname = os.path.join(odir, f)
            fhere = h5py.File(fname, 'r')
            data.append([fhere[dataset][:],fhere.attrs['TIME'],fhere.attrs['DT']])
            # xaxis = fhere['AXIS/AXIS1'][:]
            # xaxiselems = [xaxis[0]+n*data[0][2][0] for n in range(len(data[0][0]))]

        plotdata = np.zeros( (len(data),len(data[0][0])) )
        taxis = []
        for i in range(len(data)):
            plotdata[i,:] = data[i][0]
            taxis.append

    plt.figure(figsize=(8, 5))
    plt.imshow(plotdata,
               origin='lower',
               aspect='auto',
              extent=[xaxismin, xaxismax, taxismin, taxismax],
              cmap="nipy_spectral")
    if(xlim != [-1,-1]):
        plt.xlim(xlim)
    if(ylim != [-1,-1]):
        plt.ylim(ylim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)

    plt.colorbar(orientation='vertical')
    plt.title('Time vs Space')
    #plt.xlabel('$x_1 [c/\omega_p]$')
    #plt.ylabel('$p_1 [m_ec]$')

    plt.show()

    return plotdata


def wkplot(rundir='',dataset='e3',klim=[-1,-1],wlim=[-1,-1],zlim=[-1,-1],
    plotdata=[]):

    workdir = os.getcwd()
    workdir = os.path.join(workdir, rundir)

    odir = os.path.join(workdir, 'MS', 'FLD', dataset)
    files = os.listdir(odir)

    fhere = h5py.File(os.path.join(odir, files[0]), 'r')
    xaxismin = fhere['AXIS']['AXIS1'][0]
    xaxismax = fhere['AXIS']['AXIS1'][1]
    taxismin = fhere.attrs['TIME']

    fhere = h5py.File(os.path.join(odir, files[len(files)-1]), 'r')
    taxismax = fhere.attrs['TIME']

    nx = len(fhere[dataset][:])
    dx = (xaxismax-xaxismin)/nx
    nt = len(files)

    # fhere = h5py.File(odir+files[len(files)-1], 'r')
    fhere2 = h5py.File(os.path.join(odir, files[len(files)-2]), 'r')
    dt = fhere.attrs['TIME']-fhere2.attrs['TIME']
    # print(dt)

    kaxis = np.fft.fftfreq(nx, d=dx) * 2*np.pi
    waxis = np.fft.fftfreq(nt, d=dt) * 2*np.pi

    if(plotdata == []):

        data = []
        for f in files:
            fname = os.path.join(odir, f)
            fhere = h5py.File(fname, 'r')
            data.append([fhere[dataset][:],fhere.attrs['TIME'],fhere.attrs['DT']])
            # xaxis = fhere['AXIS/AXIS1'][:]
            # xaxiselems = [xaxis[0]+n*data[0][2][0] for n in range(len(data[0][0]))]

        a = np.zeros( (len(data),len(data[0][0])) )
        taxis = []
        for i in range(len(data)):
            a[i,:] = data[i][0]
            taxis.append

        plotdata = np.fliplr(np.fft.fftshift(np.fft.fft2(a)))
        plotdata = np.log(np.abs(plotdata.real))


    fig = plt.figure(figsize=(8, 8))
    plt.imshow(plotdata,
              origin='lower',
               aspect='auto',
              extent=[ min(kaxis), max(kaxis),
                     min(waxis), max(waxis) ],
              cmap="nipy_spectral")

    plt.colorbar(orientation='vertical')

    if(klim != [-1,-1]):
        plt.xlim(klim)
    if(wlim != [-1,-1]):
        plt.ylim(wlim)
    if(zlim != [-1,-1]):
        plt.clim(zlim)

    #fig.show()

    return plotdata

# def modfig(oldfig,klim=[-1,-1],wlim=[-1,-1],zlim=[-1,-1],):
#     #fig = oldfig #plt.figure(figsize=(8, 8))
#     #ax = oldfig.add_subplot(111)
#     ax = oldfig.get_axes()
#
#     if(klim != [-1,-1]):
#         ax.xlim(self,klim)
#     if(wlim != [-1,-1]):
#         ax.ylim(wlim)
#     if(zlim != [-1,-1]):
#         ax.clim(vmin=zlim[0],vmax=zlim[1])
#
#     #oldfig
#     IPython.display.display(oldfig)
#     return oldfig
