def lens1memnamegen(nmems):
    memstr=[]
    for imem in range(0,nmems,1):
        memstr1 = str(imem+1).zfill(3)
        if (imem > 34):
            memstr1 = str(imem-35+101).zfill(3)
        memstr.append(memstr1)

    return memstr

def lens2memnamegen(nmems):
    """Generate the member names for LENS2 simulations
    Input: 
      nmems = number of members
    Output:
      memstr(nmems) = an array containing nmems strings corresponding to the member names
    """

    memstr=[]
    for imem in range(0,nmems,1):
   
        if (imem < 10):
            memstr1=str(1000+imem*20+1)
            memstr2=str(imem+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)
       
        if ((imem >= 10) and (imem < 20)):
            memstr1=str(1000+(imem-10)*20+10+1)
            memstr2=str(imem+1-10).zfill(3)
            memstr.append(memstr1+'.'+memstr2)
 
        if ((imem >= 20) and (imem < 40)):
            memstr1=str(1231)
            memstr2=str(imem-20+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 40) and (imem < 60)):
            memstr1=str(1251)
            memstr2=str(imem-40+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 60) and (imem < 80)):
            memstr1=str(1281)
            memstr2=str(imem-60+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if (imem >= 80):
            memstr1=str(1301)
            memstr2=str(imem-80+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

 
    return memstr

def lens2memnamegen_first50(nmems):
    """Generate the member names for members 1-49 of LENS2
    Input: 
      nmems = number of members
    Output:
      memstr(nmems) = an array containing nmems strings corresponding to the member names
    """

    memstr=[]
    for imem in range(0,nmems,1):
   
        if (imem < 10):
            memstr1=str(1000+imem*20+1)
            memstr2=str(imem+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)
       
        if ((imem >= 10) and (imem < 20)):
            memstr1=str(1231)
            memstr2=str(imem-10+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 20) and (imem < 30)):
            memstr1=str(1251)
            memstr2=str(imem-20+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 30) and (imem < 40)):
            memstr1=str(1281)
            memstr2=str(imem-30+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if (imem >= 40):
            memstr1=str(1301)
            memstr2=str(imem-40+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

 
    return memstr





def lens2memnamegen_second50(nmems):
    """Generate the member names for members 50-100 of LENS2
    Input: 
      nmems = number of members
    Output:
      memstr(nmems) = an array containing nmems strings corresponding to the member names
    """

    memstr=[]
    for imem in range(0,nmems,1):
   
        if (imem < 10):
            memstr1=str(1000+imem*20+10+1)
            memstr2=str(imem+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)
       
        if ((imem >= 10) and (imem < 20)):
            memstr1=str(1231)
            memstr2=str(imem-10+10+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 20) and (imem < 30)):
            memstr1=str(1251)
            memstr2=str(imem-20+10+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 30) and (imem < 40)):
            memstr1=str(1281)
            memstr2=str(imem-30+10+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if (imem >= 40):
            memstr1=str(1301)
            memstr2=str(imem-40+10+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

 
    return memstr






def lens2memnamegen_temp(nmems):
    """Generate the member names for LENS2 simulations
    Input: 
      nmems = number of members
    Output:
      memstr(nmems) = an array containing nmems strings corresponding to the member names
    """

    memstr=[]
    for imem in range(0,nmems,1):
   
        if (imem < 10):
            memstr1=str(1000+imem*20+1)
            memstr2=str(imem+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)
        
        if ((imem >= 10) and (imem < 30)):
            memstr1=str(1231)
            memstr2=str(imem-10+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 30) and (imem < 50)):
            memstr1=str(1251)
            memstr2=str(imem-30+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 50) and (imem < 70)):
            memstr1=str(1281)
            memstr2=str(imem-50+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 70) and (imem < 90)):
            memstr1=str(1301)
            memstr2=str(imem-70+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

 
    return memstr





def lens2memnamegen_first50(nmems):
    """Generate the member names for LENS2 simulations
    Input: 
      nmems = number of members
    Output:
      memstr(nmems) = an array containing nmems strings corresponding to the member names
    """

    memstr=[]
    for imem in range(0,nmems,1):
   
        if (imem < 10):
            memstr1=str(1000+imem*20+1)
            memstr2=str(imem+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)
        
        if ((imem >= 10) and (imem < 20)):
            memstr1=str(1231)
            memstr2=str(imem-10+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 20) and (imem < 30)):
            memstr1=str(1251)
            memstr2=str(imem-20+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 30) and (imem < 40)):
            memstr1=str(1281)
            memstr2=str(imem-30+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

        if ((imem >= 40) and (imem < 50)):
            memstr1=str(1301)
            memstr2=str(imem-40+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

 
    return memstr





def lens2groups():
    """Generate the member names for LENS 2 simulations an output them in groups
    that correspond to the initialization type"""

    memmacro=[]
    for imem in range(0,10,1):
        memstr1=str(1000+imem*20+1)
        memstr2=str(imem+1).zfill(3)
        memmacro.append(memstr1+'.'+memstr2)
    for imem in range(0,10,1):
        memstr1=str(1000+10+(imem*20)+1)
        memstr2=str(imem+1).zfill(3)
        memmacro.append(memstr1+'.'+memstr2)

    memmicro1=[]
    for imem in range(0,20,1):
        memstr1=str(1231)
        memstr2=str(imem+1).zfill(3)
        memmicro1.append(memstr1+'.'+memstr2)

    memmicro2=[]
    for imem in range(0,20,1):
        memstr1=str(1251)
        memstr2=str(imem+1).zfill(3)
        memmicro2.append(memstr1+'.'+memstr2)

    memmicro3=[]
    for imem in range(0,20,1):
        memstr1=str(1281)
        memstr2=str(imem+1).zfill(3)
        memmicro3.append(memstr1+'.'+memstr2)

    memmicro4=[]
    for imem in range(0,20,1):
        memstr1=str(1301)
        memstr2=str(imem+1).zfill(3)
        memmicro4.append(memstr1+'.'+memstr2)

    return memmacro,memmicro1,memmicro2,memmicro3,memmicro4
























