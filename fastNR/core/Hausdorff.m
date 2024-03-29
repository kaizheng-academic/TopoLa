function [h] = Hausdorff(p1,p2)

d = zeros(length(p1),1);
for i=1:length(p1)
    t = p2-p1(i,:);
    d(i) = min(sqrt(t(:,1).^2+t(:,2).^2));    
end
hab=max(d)

for i=1:length(p2)
    t = p1-p2(i,:);
    d(i) = min(sqrt(t(:,1).^2+t(:,2).^2));    
end
hba=max(d)

h=max([hab hba])

end
