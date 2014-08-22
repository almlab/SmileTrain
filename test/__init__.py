import StringIO

def fake_fh(content=''):
    '''make a fake filehandle with this content'''
    if isinstance(content, list):
        content = "\n".join(content)
    
    fh = StringIO.StringIO()
    fh.write(content)
    fh.seek(0)
    return fh