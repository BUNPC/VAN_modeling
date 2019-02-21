function volg = getNodeVol( Vm, nodePos )

nng = size(nodePos,1);
sizeVm = size(Vm);
lstV = find(Vm>0);
Vn = zeros(sizeVm);

nodePostmp = round(nodePos);
for ii = 1:size(nodePos,1)
    % If I want branching nodes to have zero vol then exlude them here
    Vn(nodePostmp(ii,2),nodePostmp(ii,1),nodePostmp(ii,3)) = ii;
end
B = zeros(3,3,3);
B(2,2,:) = 1;
B(2,:,2) = 1;
B(:,2,2) = 1;
lst = find(Vn(lstV)==0);
nlst = length(lst);
nlst0 = nlst+10;
niter = 0;
while (nlst0-nlst)~=0
    niter = niter + 1;
    disp( sprintf('iter %d : %d voxels left',niter,nlst) )
    drawnow
    Vn2 = imdilate(Vn,B);
    Vn(lstV(lst)) = Vn2(lstV(lst));
    lst = find(Vn(lstV)==0);
    nlst0 = nlst;
    nlst = length(lst);
end  
    
% h = 10;
% hWait = waitbar( 0, 'Calculating Node Volumes (part 1 of 2)...');
% nV = length(lstV);
% for iV = 1:nV
%     waitbar(iV/nV,hWait);
%     
%     [rv(2),rv(1),rv(3)] = ind2sub(sizeVm,lstV(iV));
%     lstN = find( nodePos(:,1)-h<=rv(1) & nodePos(:,1)+h>=rv(1) & ...
%         nodePos(:,2)-h<=rv(2) & nodePos(:,2)+h>=rv(2) & ...
%         nodePos(:,3)-h<=rv(3) & nodePos(:,3)+h>=rv(3) );
%     if ~isempty(lstN)
%         [foo,iN] = min( (sum(((ones(length(lstN),1)*rv)-nodePos(lstN,:)).^2,2) ).^0.5 );
%         Vn(lstV(iV)) = lstN(iN);
%     end
% %    [foo, iN] = min( (sum(((ones(nng,1)*rv)-nodePos).^2,2) ).^0.5 );
% %    Vn(lstV(iV)) = iN;
% end
% close(hWait);   



hWait = waitbar( 0, 'Calculating Node Volumes (part 2 of 2)...');
volg = zeros(nng,1);
for ii=1:nng
    waitbar(ii/nng,hWait);
%    volg(ii) = length(find(Vn(lstV)==ii));
    lst = find(Vn(lstV)==ii);
    volg(ii) = length(lst);
    lstV(lst) = [];
end
close(hWait);   


