
from uuid import uuid4

class AnalysisHandler(object):

    def __init__(self) -> None:
        self.job_id = uuid4()
