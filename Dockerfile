# app/Dockerfile

FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

#RUN git clone https://github.com/streamlit/streamlit-example.git .

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

#HEALTHCHECK CMD curl --fail http://localhost:8502/_stcore/health

#ENTRYPOINT ["streamlit", "run", "main.py", "--server.port=8502", "--server.address=0.0.0.0"]

CMD ["watchmedo", "auto-restart", "--directory", ".", "--pattern", "*.py", "--recursive", "--", "streamlit", "run", "main.py", "--server.port", "8502", "--server.address", "0.0.0.0"]
