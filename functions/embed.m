function Pic_out=embed(Pic_in, target_size, padval)
%Pic_out=embed(Pic_in, target_size, padval)
%embeds a vector, array or cube inside a bigger canvas the size of which is
%determined by target_size (target size can be a vector)
%padval...padding value


% %test parameters
% padval=0;
% %2D
% Pic_in=zeros(26,26); Pic_in(14,14)=1; target_size=[13 13]; %2D - even to odd, getting smaller
% % Pic_in=zeros(26,26); Pic_in(14,14)=1; target_size=[14 14]; %2D - even to even, getting smaller
% % Pic_in=zeros(25,25); Pic_in(13,13)=1; target_size=[12 12]; %2D - odd to even, getting smaller
% Pic_in=zeros(25,25); Pic_in(13,13)=1; target_size=[13 13]; %2D - odd to odd, getting smaller
% 
% %Pic_in=zeros(26,26); Pic_in(14,14)=1; target_size=[31 31]; %2D - even to odd, getting larger
% %Pic_in=zeros(26,26); Pic_in(14,14)=1; target_size=[34 34]; %2D - even to even, getting larger
% %Pic_in=zeros(25,25); Pic_in(13,13)=1; target_size=[32 32]; %2D - odd to even, getting larger
% %Pic_in=zeros(25,25); Pic_in(13,13)=1; target_size=[33 33]; %2D - odd to odd, getting larger




if ndims(Pic_in)==2 && min(size(Pic_in))==1; %1D
    
    n_1=ceil((target_size-length(Pic_in))/2);
    n_2=floor((target_size-length(Pic_in))/2);
    
    if size(Pic_in,1)==1; %if row-vector
        tmp=padarray(Pic_in,[0 n_1],padval,'pre');
        Pic_out=padarray(tmp,[0 n_2],padval,'post');
    else %if column-vector
        tmp=padarray(Pic_in,[n_1 0],padval,'pre');
        Pic_out=padarray(tmp,[n_2 0],padval,'post');
    end
    
    
%-------------------------------------------------------------------------

        
elseif ndims(Pic_in)==3; %3D
    
    if length(target_size)==1; target_size=[target_size target_size target_size]; end;

    siz=size(Pic_in);
    size_diff=target_size(1:3)-siz(1:3);
    
     %adapt number of rows:
        if size_diff(1)>0; %if height of Pic_in is smaller than the target_size
            
            if mod(target_size(1),1)==0;
                tmp=padarray(Pic_in,[floor(size_diff(1)/2) 0 0],padval,'post');
                Pic_out2=padarray(tmp,[ceil(size_diff(1)/2) 0 0],padval,'pre');
            else
                tmp=padarray(Pic_in,[ceil(size_diff(1)/2) 0 0],padval,'post');
                Pic_out2=padarray(tmp,[floor(size_diff(1)/2) 0 0],padval,'pre');
            end
            
        else
            if mod(size(Pic_in,1),2)==0;
                vec=(ceil(abs(size_diff(1))/2)+1):ceil(abs(size_diff(1))/2)+target_size(1);
            else
                vec=(floor(abs(size_diff(1))/2)+1):floor(abs(size_diff(1))/2)+target_size(1);
            end
            Pic_out2=Pic_in(vec,:,:);
        end

    %adapt number of columns:
        if size_diff(2)>0; %if height of Pic_in is smaller than the target_size
            
            if mod(target_size(2),1)==0;
                tmp=padarray(Pic_out2,[0 floor(size_diff(2)/2) 0],padval,'post');
                Pic_out1=padarray(tmp,[0 ceil(size_diff(2)/2) 0],padval,'pre');
            else
                tmp=padarray(Pic_out2,[0 ceil(size_diff(2)/2) 0],padval,'post');
                Pic_out1=padarray(tmp,[0 floor(size_diff(2)/2) 0],padval,'pre');
            end
            
        else
            if mod(size(Pic_in,2),2)==0;
                vec=(ceil(abs(size_diff(2))/2)+1):ceil(abs(size_diff(2))/2)+target_size(2);
            else
                vec=(floor(abs(size_diff(2))/2)+1):floor(abs(size_diff(2))/2)+target_size(2);
            end
            Pic_out1=Pic_out2(:,vec,:);
        end
    
    %adapt number of pages:
        if size_diff(3)>0; %if height of Pic_in is smaller than the target_size
            
            if mod(target_size(2),1)==0;
                tmp=padarray(Pic_out1,[0 0 floor(size_diff(3)/2)],padval,'post');
                Pic_out=padarray(tmp,[0 0 ceil(size_diff(3)/2)],padval,'pre');
            else
                tmp=padarray(Pic_out1,[0 0 ceil(size_diff(3)/2)],padval,'post');
                Pic_out=padarray(tmp,[0 0 floor(size_diff(3)/2)],padval,'pre');
            end
            
        else
            if mod(size(Pic_in,3),2)==0;
                vec=(ceil(abs(size_diff(3))/2)+1):ceil(abs(size_diff(3))/2)+target_size(3);
            else
                vec=(floor(abs(size_diff(3))/2)+1):floor(abs(size_diff(3))/2)+target_size(3);
            end
            Pic_out=Pic_out1(:,:,vec);
        end
        
%-------------------------------------------------------------------------        
        
else %2D
        
    %pads images to a defined size; the frame around is filled with padval
    %targetsize can be a vector: [rows colmns]

    if length(target_size)==1; target_size=[target_size target_size]; end;
    %if length(size(Pic_in))>=3; Pic_in=sum(Pic_in,3)/size(Pic_in,3); end;

    siz=size(Pic_in);
    size_diff=target_size(1:2)-siz(1:2);

    %adapt number of rows:
        if size_diff(1)>0 %if height of Pic_in is smaller than the target_size

            if mod(target_size(1),2)==0
                tmp=padarray(Pic_in,[floor(size_diff(1)/2) 0],padval,'post');
                Pic_out1=padarray(tmp,[ceil(size_diff(1)/2) 0],padval,'pre');
            else
                tmp=padarray(Pic_in,[ceil(size_diff(1)/2) 0],padval,'post');
                Pic_out1=padarray(tmp,[floor(size_diff(1)/2) 0],padval,'pre');
            end
            
        else %if height of input is larger than that of output
            if mod(size(Pic_in,1),2)==0                 
                vec=(ceil(abs(size_diff(1))/2)+1):(ceil(abs(size_diff(1))/2))+target_size(1);
            else
                vec=(floor(abs(size_diff(1))/2)+1):(floor(abs(size_diff(1))/2))+target_size(1);
            end
            Pic_out1=Pic_in(vec,:);
        end


    %adapt number of columns:
        if size_diff(2)>0 %if height of Pic_in is smaller than the target_size

            if mod(target_size(2),2)==0
                tmp=padarray(Pic_out1,[0 floor(size_diff(2)/2)],padval,'post');
                Pic_out=padarray(tmp,[0 ceil(size_diff(2)/2)],padval,'pre');
            else
                tmp=padarray(Pic_out1,[0 ceil(size_diff(2)/2)],padval,'post');
                Pic_out=padarray(tmp,[0 floor(size_diff(2)/2)],padval,'pre');
            end
            

        else
            if mod(size(Pic_in,2),2)==0
                vec=(ceil(abs(size_diff(2))/2)+1):(ceil(abs(size_diff(2))/2))+target_size(2);
            else
                vec=(floor(abs(size_diff(2))/2)+1):(floor(abs(size_diff(2))/2))+target_size(2);
            end
            Pic_out=Pic_out1(:,vec);
        end
end

%imagesc(Pic_out); %only for tests
    
    