function ptha=wdijkstra(C,w,nc)
% find pathes with dijkstra's algotithm  http://en.wikipedia.org/wiki/Dijkstra's_algorithm
% C -connection matrix
% w- wights matrix
% nc- for wich node find shortest paths
n=size(w,1); % number of nodes

ptha=cell(n,1); % paths, empty at beginig
for nct=1:n
    ptha{nct}=nc; % add start from main node
end
%ptha{nc} - path to itself will keep consisted of nc-node only

% markers, initial values:
m=Inf(1,n);
m(nc)=0;


used=false(1,n); % all vertecies not used at begining

while true
    if all(used)
        break; % if no more unused nodes => finish
    end
    nu=find(~used);
    [mmn0 u0]=min(m(nu)); % minimal makrer node is next
    u=nu(u0); % not used node with minimal marker is current node
    f=find(C(u,:)&(~used)); % connected and not used nighbours 

    
    % pathes:
    for fc=1:length(f) % for each reachible neighbour
        f1=f(fc);
        ov1=m(f1); % old value
        nv1=m(u)+w(u,f1); % new valu
        if nv1<ov1
            ptha{f1}=[ptha{u} f1]; % replce old path
            m(f1)=m(u)+w(u,f1);
        end
        
    end
    used(u)=true; % all neighbours processed => current node become used
    
    
    
end