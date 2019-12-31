FROM python:3.7.6-buster

COPY requirements.txt /tmp/

RUN apt-get update && apt-get -y install \
	gcc 

RUN pip3 install -r /tmp/requirements.txt

EXPOSE 8888

CMD ["jupyter","notebook","--allow-root","--ip","0.0.0.0"]
