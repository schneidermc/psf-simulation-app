function stack_out = remove_BG(stack)
    [Nx, Ny, Nz ] = size(stack);
    BG(1,1,:)=min([squeeze(mean(stack(1,:,:),2)), squeeze(mean(stack(end,:,:),2)), squeeze(mean(stack(:,1,:),1)), squeeze(mean(stack(:,end,:),1))],[],2);
    stack_out = stack - repmat(BG,[Nx Nx 1]); %clear I2;
    stack_out(stack_out < 0) = 0; %ensuring positivity    