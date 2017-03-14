import xml.etree.ElementTree as ET
import os.path as path

class Param(object):
    """A python way to implement a Fortran namelist.

    Parameters are read from a xml file

    The xml file contains defaults value and list of available values
    for parameters who take only a limited set of values
    """
    def __init__(self,xmlfile='default.xml',verbose=True):
        import grid
        d = path.dirname(grid.__file__)
        #print('fluid2d main core is in %s'%d)
        if xmlfile=='default.xml':
            self.xmlfile = d+'/default.xml'
        else:
            self.xmlfile = xmlfile
        self.verbose = verbose
        self.set_attr_from_xml()
        
    def __str__(self):
        print('Parameters are:')
        for p in self.list_param:
            print(' - %s = %s'%(p,str(getattr(self,p))))
        print('default parameters come from: %s'%self.xmlfile)
        return ''
            
    def set_attr_from_xml(self):
        """ read the xml and assign default value to the attributes """
        tree = ET.parse(self.xmlfile)
        root = tree.getroot()
        self.list_param=[]
        for p in root.iter('name'):
            name=p.attrib['value']
            typ = p.find('type').text.split(',')[0]    
            tex = p.find('default').text
            val = eval("%s('%s')"%(typ,tex))    
            if(typ=='bool'):
                val = tex in ['True','true']
            setattr(self,name,val)
            self.list_param.append(name)

    def check(self):
        """ check whether the values are within the available values """
        tree = ET.parse(self.xmlfile)
        root = tree.getroot()
        generalok = True
        for p in root.iter('name'):
            name=p.attrib['value']
            typ = p.find('type').text.split(',')[0]    
            avail = p.find('avail')
            if hasattr(avail,'text'):
                list_avail = avail.text.replace(' ','').split(',')
                val = getattr(self,name)
                ok = (val in list_avail)
                if ok:
                    if self.verbose:
                        print('value of %s is ok'%name)
                else:
                    generalok = False
                    if self.verbose:
                        print('value of %s is wrong, should be one of [%s]'%(name,avail.text))
        return generalok                
        
    def copy(self,obj,list_param):
        """ copy attributes listed in list_param to obj
        
        On output it returns missing attributes
        """
        missing=[]
        for k in list_param:
            if hasattr(self,k):
                setattr(obj,k,getattr(self,k))
            else:
                missing.append(k)
        return missing
