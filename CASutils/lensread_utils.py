

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
            memstr1=str(1401)
            memstr2=str(imem-40+1).zfill(3)
            memstr.append(memstr1+'.'+memstr2)

 
    return memstr
