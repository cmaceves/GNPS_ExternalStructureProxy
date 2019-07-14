build:
	docker build -t externalstructureproxy .

server:
	docker run --rm -d -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy /app/run_server.sh

interactive:
	docker run --rm -it -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy /app/run_server.sh

dev-server:
	docker run --rm -it -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy /app/run_dev_server.sh

bash:
	docker run --rm -it -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy bash