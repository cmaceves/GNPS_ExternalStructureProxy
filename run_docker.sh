docker rm npatlasgnpsproxy
#docker run -d  -p 5010:5010 --name npatlasgnpsproxy gnpsmetadata /app/run_production_server.sh
#docker run -d -p 5010:5010 --name npatlasgnpsproxy npatlasgnpsproxy /app/run_server.sh
docker run -it -p 5010:5010 --name npatlasgnpsproxy npatlasgnpsproxy bash
