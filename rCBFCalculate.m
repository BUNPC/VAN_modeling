function rCBF=rCBFCalculate(savefolder,VesselDiPercent) % to calculate rCBF. By Xiaojun Cheng July 2017.

load([savefolder,'mesh.mat']);
load([savefolder,'dilate_vessel_NCES_Sigmoid',num2str(VesselDiPercent),'.mat']);
nn=1;
FlowInNodest = find(im2.Adv.nodeIn==1); % find the artery
for ii=1:length(FlowInNodest)
if im2.nodeType(FlowInNodest(ii))==1 %artery
        FlowInNodes(nn)=FlowInNodest(ii);nn=nn+1;
    
end
end


L=length(FlowInNodes);

for ii=1:L
    edgeindex=find(im2.nodeEdges(:,1)==FlowInNodes(ii));
    if length(edgeindex)==0
        edgeindex=find(im2.nodeEdges(:,2)==FlowInNodes(ii));
    end
    if edgeindex>0
        CBF_base(ii)=Fedges(edgeindex,1);
        CBF_final(ii)=Fedges(edgeindex,end);
    else
        CBF_base(ii)=0;
        CBF_final(ii)=0;
        
    end
    
end
rCBF=sum(abs(CBF_final))/sum(abs(CBF_base));

end