class FakeFile:
    def __init__(self, lines):
        self.lines = lines
        
    def readline(self):
        if len(self.lines) > 0:
            return self.lines.pop(0) + "\n"
        else:
            return ''