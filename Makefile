build:
	docker build -t externalstructureproxy .

server:
	docker run --rm -d -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy /app/run_production_server.sh

interactive:
	docker run --rm -it -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy /app/run_production_server.sh

dev-server:
	docker run --rm -it -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy /app/run_server.sh

bash:
	docker run --rm -it -p 5010:5000 --rm --name externalstructureproxy externalstructureproxy bash


server-compose:
	docker-compose build
	docker-compose up -d

server-compose-interactive:
	docker-compose build
	docker-compose up

server-compose-production:
	docker-compose build
	docker-compose -f docker-compose.yml -f docker-compose-production.yml up -d
