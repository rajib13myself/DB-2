clear
ls
mkdir -p /home/ubuntu/cluster-keys
ls
pwd
whoami
cd ssh
cd .ssh
ls
exit
celar
clea
clear
ls
ssh-keygen -t rsa
rm /home/ubuntu/.ssh/id_rsa
rm /home/ubuntu/.ssh/id_rsa.pub
ls /home/ubuntu/.ssh
ssh-keygen -t rsa
ls
cd cluster-keys/
ls
cd ..
cd .ssh
ls
cd ..
pwd
sudo cat /home/ubuntu/.ssh/cluster-key.pub
hostname -I
sudo cat /var/lib/rancher/k3s/server/node-token
clear
sudo systemctl status k3s
curl -sfL https://get.k3s.io | sh -
sudo systemctl status k3s
sudo cat /var/lib/rancher/k3s/server/node-token
ls
hostname -I
ls
mv 'UPPMAX 2024_1-4-openrc.sh' UPPMAX\ 2024_1-4-openrc.sh
mv 'UPPMAX 2024_1-4-openrc.sh' UPPMAX\ 2024_1-4-openrc1.sh
mv UPPMAX 2024_1-4-openrc1.sh UPPMAX\ 2024_1-4-openrc.sh
mv 'UPPMAX 2024_1-4-openrc1.sh' UPPMAX\ 2024_1-4-openrc.sh
ls
clear
source UPPMAX\ 2024_1-4-openrc.sh
ls
sudo apt update
sudo apt install python3-pip
python3 start_instance.py 
ls
exit
clear
source UPPMAX\ 2024_1-4-openrc.sh
python3 start_instance.py 
pip3 install python-novaclient
sudo apt-get update
sudo apt-get install tree
python3 start_instance.py 
pip3 install python-keystoneclient
python3 start_instance.py
nano start_instance.py
ls
nano start_instance.py
python3 start_instance.py
nano start_instance.py
python3 start_instance.py
HOSTNAME i
HOSTNAME I
Hostname I
hostname -I
ping 192.168.2.79
ping 192.168.2.48
clear
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.79
ls
source UPPMAX\ 2024_1-4-openrc.sh
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.79
exit
ls
clear
sudo apt-get update
sudo apt-get install tree
tree
clear
sudo apt install -y python3-pip
source UPPMAX\ 2024_1-4-openrc.sh
ls
source UPPMAX\ 2024_1-4-openrc.sh
ls
pip3 install python-novaclient
pip3 install python-keystoneclient
ls
cat start_instance.py 
nano start_instance.py 
python3 start_instance.py 
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.79
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.68
ping 192.168.2.24
ping 192.168.2.68
exit
ls
ping 192.168.2.24
ping 192.168.2.68
ls
cd .ssh
cp /home/ubuntu/RD_SSH.pem .
ls
cd ..
ssh -i .ssh/RD_SSH.pem ubuntu@192.168.2.68
cd .ssh
pwd
cd ..
ssh -i /home/ubuntu/.ssh/RD_SSH.pem ubuntu@192.168.2.68
chmod 400 ~/.ssh/RD_SSH.pem
ssh-add ~/.ssh/RD_SSH.pem
scp -i ~/.ssh/RD_SSH.pem RD_SSH.pem ubuntu@192.168.2.68:/home/ubuntu/.ssh/
scp -i ~/.ssh/RD_SSH.pem /home/ubuntu/RD_SSH.pem ubuntu@192.168.2.68:/home/ubuntu/.ssh/
exit
ls
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
ssh-keygen -f "/home/ubuntu/.ssh/known_hosts" -R "192.168.2.24"
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
cd .ssh/
ls
cd ..
scp -i /home/ubuntu/.ssh/cluster-key -r ./K8-Task5/ appuser@192.168.2.24
scp -i /home/ubuntu/.ssh/cluster-key -r ./K8-Task5/ appuser@192.168.2.24:.
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
exit]
exit
ls
scp -i ~/.ssh/RD_SSH.pem RD_SSH.pem ubuntu@192.168.2.68:/home/ubuntu/.ssh/
ping 192.168.2.68
scp -i ~/.ssh/RD_SSH.pem /home/ubuntu/RD_SSH.pem ubuntu@192.168.2.68:/home/ubuntu/.ssh/
source UPPMAX\ 2024_1-4-openrc.sh
scp -i ~/.ssh/RD_SSH.pem /home/ubuntu/RD_SSH.pem ubuntu@192.168.2.68:/home/ubuntu/.ssh/
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
nano cloud-config-master.yml 
ssh -i /home/ubuntu/.ssh/cluster-key appuser@<Worker-VM-IP>
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
sudo cat /home/ubuntu/.ssh/cluster-key.pub
hostname -I
sudo cat /var/lib/rancher/k3s/server/node-token
clear
ssh -i /home/ubuntu/.ssh/cluster-key appuser@<Worker-VM-IP>
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.68
cat /var/lib/cloud/instance/scripts/runcmd
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.68
clear
ipconfig
sudo apt install net-tools
ipconfig
ifconfig | grep "inet "
sudo k3s kubectl get token --print-join-command
sudo cat /var/lib/rancher/k3s/server/node-token
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.68
ssh-keygen -f "/home/ubuntu/.ssh/known_hosts" -R "192.168.2.68"
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.68
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
ls
cd K8-Task5/
l
cd ..
scp -i .ssh/RD_SSH.pem -r ./K8-Task5/ ubuntu@192.168.2.68:/home/ubuntu/
scp -i /home/ubuntu/.ssh/cluster-key -r ./K8-Task5/ appuser@192.168.2.68:/home/ubuntu/
ls
docker build -t rajib43datta/test-python-app:latest .
sudo apt  install docker.io
ls
sudo docker --version sion
docker build -t rajib43datta/test-python-app:latest .
nano Dockerfile
nano app.py
nano deployment.yaml
nano requirement.txt
mv requirement.txt requirements.txt 
nano service.yaml
ls
docker build -t rajib43datta/test-python-app:latest .
sudo docker build -t rajib43datta/test-python-app:latest .
docker login
clear
docker login
sudo docker login
docker push your-dockerhub-username/your-python-app:latest
sudo docker push your-dockerhub-username/your-python-app:latest
sudo docker push rajib43datta/your-python-app:latest
docker push rajib43datta/test-python-app:latest
sudo docker push rajib43datta/test-python-app:latest
nano namespace.yaml
ls
kubectl apply -f namespace.yaml
sudo kubectl apply -f namespace.yaml
sudo kubectl apply -f deployment.yaml
sudo kubectl apply -f service.yaml
sudo kubectl get all
sudo kubectl logs web-deployment-f87ddb556-dt895
sudo kubectl log web-deployment
kubectl log web-deployment
kubectl log web-deployment top
kubectl log top
kubectl get pods -n my-namespace
kubectl get svc -n my-namespace
kubectl get svc -n datta-namespace
sudo kubectl get svc -n datta-namespace
sudo kubectl get pods -n datta-namespace
./k3s-example
mkdir ks3-example
ls
mv namespace.yaml ./ks3-example/
mv deployment.yaml requirements.txt service.yaml ./ks3-example/
ls
sudo kubectl apply -f ./k3s-example/namespace.yaml
sudo kubectl apply -f ~/k3s-example/namespace.yaml
sudo kubectl apply -f /home/ubuntu/k3s-example/namespace.yaml
ls
sudo kubectl apply -f /home/ubuntu/ks3-example/namespace.yaml
cd ks3-example/
ls
cat requirements.txt 
cat deployment.yaml 
cd ..
sudo kubectl apply -f /home/ubuntu/ks3-example/deployment.yaml
nano ./ks3-example/deployment.yaml 
nano ./ks3-example/service.yaml 
sudo kubectl apply -f /home/ubuntu/ks3-example/deployment.yaml
sudo kubectl apply -f /home/ubuntu/ks3-example/service.yaml
sudo kubectl get pods -n datta-namespace
kubectl get svc -n datta-namespace
sudo kubectl get svc -n datta-namespace
kubectl logs <pod-name> -n datta-namespace
kubectl logs web-deployment-f87ddb556-jxm2c -n datta-namespace
sudo kubectl logs web-deployment-f87ddb556-jxm2c -n datta-namespace
nano ./ks3-example/service.yaml 
nano ./ks3-example/deployment.yaml 
nano Dockerfile 
sudo kubectl apply -f ./k3s-example/
sudo kubectl apply -f ./kss-example/
sudo kubectl apply -f /home/ubuntu/kss-example/
ls
sudo kubectl apply -f /home/ubuntu/ks3-example/
nano ./ks3-example/deployment.yaml 
sudo kubectl apply -f /home/ubuntu/kss-example/
sudo kubectl apply -f /home/ubuntu/ks3-example/
sudo kubectl get pods -n datta-namespace
kubectl get svc -n datta-namespace
sudo kubectl get svc -n datta-namespace
kubectl logs <pod-name> -n datta-namespace
sudo kubectl logs web-deployment-c7c79ddb-jzskm -n datta-namespace
sudo kubectl describe web-deployment-c7c79ddb-jzskm -n datta-namespace
sudo kubectl get svc -n datta-namespace
sudo kubectl get pods -n datta-namespace
sudo kubectl describe svc web-service -n datta-namespace
sudo kubectl get pods -n datta-namespace
sudo kubectl describe pods web-deployment-c7c79ddb-jzskm -n datta-namespace
kubectl get all
sudo kubectl get all
curl -sfL https://get.k3s.io | sh -
sudo systemctl status k3s
sudo cat /var/lib/rancher/k3s/server/node-token
hostname -I
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.68
ssh-keygen -f "/home/ubuntu/.ssh/known_hosts" -R "192.168.2.68"
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.68
sudo kubectl get nodes
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
sudo kubectl get nodes
sudo kubectl apply -f /home/ubuntu/ks3-example/
ls
sudo kubectl get nodes
kubectl get pods --all-namespaces
sudo kubectl get pods --all-namespaces
mongo --host 192.168.2.24
sudo mongo --host 192.168.2.24
sudo apt install mongodb-clients
sudo mongo --host 192.168.2.24
exir
exit
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
mongo --host 192.168.2.24
ssh -i /home/ubuntu/.ssh/cluster-key appuser@192.168.2.24
mongo --host 192.168.2.24
exit
