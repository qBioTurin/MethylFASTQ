FROM python:3.8-slim

RUN pip install --upgrade pip && \ 
    pip install --no-cache-dir biopython

COPY src/*.py /src/

ENTRYPOINT ["/src/methylFASTQ.py"]