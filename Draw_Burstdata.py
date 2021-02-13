import urllib
import urllib2
import re
import os
#dir(re)
from socket import timeout
from socket import error as SocketError
import subprocess
import time
import sys
import ssl
from functools import wraps
def sslwrap(func):
    @wraps(func)
    def bar(*args, **kw):
        kw['ssl_version'] = ssl.PROTOCOL_TLSv1
        return func(*args, **kw)
    return bar
ssl.wrap_socket = sslwrap(ssl.wrap_socket)
year=sys.argv[2:]
yourdir=sys.argv[1]
rooturl="https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts/"
if not os.path.exists(yourdir+'BurstData') :
    os.makedirs(yourdir+'BurstData')
def format_size(bytes):
    try:
        bytes = float(bytes)
        kb = bytes / 1024
    except:
        print("format is wrong")
        return "Error"
    if kb >= 1024:
        M = kb / 1024
        if M >= 1024:
            G = M / 1024
            return "%.1fG" % (G)
        else:
            return "%.1fM" % (M)
    else:
        return "%.1fK" % (kb)
def reporthook(a,b,c):
    """
        display downloading proces
        :param a: num of blocks downloaded
        :param b: size of each block
        :param c: total size
        :return: None
     """
#     speed = (a * b) / (time.time() - start_time)
# #     #speed_str = " Speed: %.2f" % format_size(speed)
#     speed_str = " %s/s" % format_size(speed)
    dowsize=(a*b)
    Downsize="%s" %format_size(dowsize)
    totalsize="%s" %format_size(c)
    #print(totalsize)
    per = 100.0 * a * b / c
    if per > 100 :
        per = 100
    #sys.stdout.write('%.2f%%' % per)
    a=per/5
    b=20-per
    sys.stdout.write('\r%s %s %s%s% d%%'
                     %(filename , Downsize,int(a)*'#',int(b)*' ',per))
    #ys.stdout.write('\r%s'%Downsize)
    sys.stdout.flush()
#     if per == 100:
#         sys.stdout.write('\n This file has been downloaded')
def reporthook1(a,b,c):
    """
        display downloading proces
        :param a: num of blocks downloaded
        :param b: size of each block
        :param c: total size
        :return: None
     """
#     speed = (a * b) / (time.time() - start_time)
# #     #speed_str = " Speed: %.2f" % format_size(speed)
#     speed_str = " %s/s" % format_size(speed)
    dowsize=(a*b)
    Downsize="%s" %format_size(dowsize)
    totalsize="%s" %format_size(c)
    #print(totalsize)
    per = 100.0 * a * b / c
    if per > 100 :
        per = 100
    #sys.stdout.write('%.2f%%' % per)
    a=per/5
    b=20-per
    sys.stdout.write('\r%s %s %s%s% d%%'
                     %(filename1 , Downsize,int(a)*'#',int(b)*' ',per))
    #ys.stdout.write('\r%s'%Downsize)
    sys.stdout.flush()
#     if per == 100:
#         sys.stdout.write('\n This file has been downloaded')
def reporthook2(a,b,c):
    """
        display downloading proces
        :param a: num of blocks downloaded
        :param b: size of each block
        :param c: total size
        :return: None
     """
#     speed = (a * b) / (time.time() - start_time)
# #     #speed_str = " Speed: %.2f" % format_size(speed)
#     speed_str = " %s/s" % format_size(speed)
    dowsize=(a*b)
    Downsize="%s" %format_size(dowsize)
    totalsize="%s" %format_size(c)
    #print(totalsize)
    per = 100.0 * a * b / c
    if per > 100 :
        per = 100
    #sys.stdout.write('%.2f%%' % per)
    a=per/5
    b=20-per
    sys.stdout.write('\r%s %s %s%s% d%%'
                     %(filename2 , Downsize,int(a)*'#',int(b)*' ',per))
    #ys.stdout.write('\r%s'%Downsize)
    sys.stdout.flush()
#     if per == 100:
#         sys.stdout.write('\n This file has been downloaded')


for year in year:    
    if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'):
        os.makedirs(yourdir+'BurstData'+'/'+year+'/')
    try:
        webcontent=urllib.urlopen(rooturl+year+"/").read()
        pat='>bn(.*?)/</a> '
        burstinfo=re.compile(pat).findall(webcontent)
        burstinfo=iter('bn'+i for i in burstinfo)
        while True:
             try:
                a=0
                b=0
                c=0
                burstdir=next(burstinfo)
                if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'+burstdir):
                    os.makedirs(yourdir+'BurstData'+'/'+year+'/'+burstdir)
                else:
                    print(yourdir+'BurstData'+'/'+year+'/'+burstdir+"exists")
                print("start downloading "+burstdir+'\n')
                foldercurrent=urllib.urlopen(rooturl+year+"/"
                    +burstdir+"/current/").read()
                pat1='>g(.*?)</a> '
                filecurrent=re.compile(pat1).findall(foldercurrent)
                #print (fileinfo)
                filecurrent=iter('g'+file for file in filecurrent)
                #print (next(fileinfo))
                #print(filename)
                folderprevious=urllib.urlopen(rooturl+year+"/"
                    +burstdir+"/previous/").read()
                fileprevious=re.compile(pat1).findall(folderprevious)
                fileprevious=iter('g'+file for file in fileprevious)
                folderquicklook=urllib.urlopen(rooturl+year+"/"
                    +burstdir+"/quicklook/").read()
                filequicklook=re.compile(pat1).findall(folderquicklook)
                filequicklook=iter('g'+file for file in filequicklook)
                while True:
                    try:
                        if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/current/'):
                            os.makedirs(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/current/')
                        else:
                            if a == 0:
                                print(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/current/'+' exists'+'\n')
                        filename=next(filecurrent)
                        filenameurl='%s%s%s%s%s%s'%(rooturl,year,"/",burstdir,"/current/",filename)
                        if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/current/'+filename):           
                            a=a+1
                            #print(a,filenameurl,year+"/"+burstdir+"/current/"+filename)
                            urllib.urlretrieve(filenameurl,yourdir+'BurstData'+'/'+
                                    year+'/'+burstdir+'/current/'+filename,reporthook) 
                        else:
                            print(filename+'exists')
                           # continue
                    except timeout:
                        print("==> Timeout")
                    except StopIteration:
                        break
                while True:
                    try:
                        if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/previous/'):
                            os.makedirs(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/previous/')
                        else:
                            if b == 0:
                                print(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/previous/'+' exists'+'\n')
                        filename1=next(fileprevious)
                        filename1url='%s%s%s%s%s%s'%(rooturl,year,"/",burstdir,"/previous/",filename1)
                        if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/previous/'+filename1):
                            b=b+1
                            #print(b,filename1url,year+"/"+burstdir+"/previous/"+filename1)
                            urllib.urlretrieve(filename1url,yourdir+'BurstData'+'/'+
                                     year+'/'+burstdir+'/previous/'+filename1,reporthook1)
                        else:
                            print(filename1+'exists')
                           # continue
                    except timeout:
                        print("==> Timeout")
                    except StopIteration:
                        break
                while True:
                    try:
                        if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/quicklook/'):
                            os.makedirs(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/quicklook/')
                        else:
                            if c == 0:
                                print(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/quicklook/'+' exists'+'\n')
                        filename2=next(filequicklook)
                        filename2url='%s%s%s%s%s%s'%(rooturl,year,"/",burstdir,"/quicklook/",filename2)
                        if not os.path.exists(yourdir+'BurstData'+'/'+year+'/'+burstdir+'/quicklook/'+filename2):
                            c=c+1
                            #print(c,filename2url,year+"/"+burstdir+"/quicklook/"+filename2)
                            urllib.urlretrieve(filename2url,yourdir+'BurstData'+'/'+
                                     year+'/'+burstdir+'/quicklook/'+filename2,reporthook2)
                        else:
                            print(filename2+'exists')
                           # continue
                    except timeout:
                        print("==> Timeout")
                    except StopIteration:
                        break     
             except timeout:
                print("==> Timeout")
             except StopIteration:
                break   
                #u=urllib.urlretrieve(filenameurl,filename,reporthook)               
    except timeout:
        print("==> Timeout")       
