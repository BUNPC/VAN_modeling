function [nB,im] = nBupdate( im )

if im.nBflag==1
    nNodes = size(im.nodePos,1);
    nB=zeros(1,nNodes);
    hwait = waitbar( 0,'Updating nB' );
    for ii=1:nNodes
        nB(ii)=length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
    end
    im.nB = nB;
    close(hwait);
else
    nB = im.nB;
end
