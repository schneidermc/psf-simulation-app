function Pic_out = embedArray(Pic_in, target_size, padval)
    % Embeds a vector, array or cube inside a bigger canvas the size of which
    % is determined by target_size (target size can be a vector)
    % padval...padding value

    if isvector(Pic_in) %1D
        n_1 = ceil((target_size-length(Pic_in))/2);
        n_2 = floor((target_size-length(Pic_in))/2);

        if size(Pic_in,1)==1 %if row-vector
            tmp = padarray(Pic_in,[0 n_1],padval,'pre');
            Pic_out = padarray(tmp,[0 n_2],padval,'post');
        else %if column-vector
            tmp = padarray(Pic_in,[n_1 0],padval,'pre');
            Pic_out = padarray(tmp,[n_2 0],padval,'post');
        end

    elseif ndims(Pic_in)==3 %3D

        if length(target_size)==1
            target_size = [target_size target_size target_size];
        end

        siz=size(Pic_in);
        size_diff=target_size-siz(1:3);

        %adapt number of rows:
        if size_diff(1)>0 %if height of Pic_in is smaller than the target_size
            if mod(target_size(1),1)==0
                tmp=padarray(Pic_in,[floor(size_diff(1)/2) 0 0],padval,'post');
                Pic_out2=padarray(tmp,[ceil(size_diff(1)/2) 0 0],padval,'pre');
            else
                tmp=padarray(Pic_in,[ceil(size_diff(1)/2) 0 0],padval,'post');
                Pic_out2=padarray(tmp,[floor(size_diff(1)/2) 0 0],padval,'pre');
            end
        else
            if mod(size(Pic_in,1),2)==0
                vec=(ceil(abs(size_diff(1))/2)+1):ceil(abs(size_diff(1))/2)+target_size(1);
            else
                vec=(floor(abs(size_diff(1))/2)+1):floor(abs(size_diff(1))/2)+target_size(1);
            end
            Pic_out2=Pic_in(vec,:,:);
        end

        %adapt number of columns:
        if size_diff(2)>0 %if height of Pic_in is smaller than the target_size
            if mod(target_size(2),1)==0
                tmp=padarray(Pic_out2,[0 floor(size_diff(2)/2) 0],padval,'post');
                Pic_out1=padarray(tmp,[0 ceil(size_diff(2)/2) 0],padval,'pre');
            else
                tmp=padarray(Pic_out2,[0 ceil(size_diff(2)/2) 0],padval,'post');
                Pic_out1=padarray(tmp,[0 floor(size_diff(2)/2) 0],padval,'pre');
            end
        else
            if mod(size(Pic_in,2),2)==0
                vec=(ceil(abs(size_diff(2))/2)+1):ceil(abs(size_diff(2))/2)+target_size(2);
            else
                vec=(floor(abs(size_diff(2))/2)+1):floor(abs(size_diff(2))/2)+target_size(2);
            end
            Pic_out1=Pic_out2(:,vec,:);
        end

        %adapt number of pages:
        if size_diff(3)>0 %if height of Pic_in is smaller than the target_size
            if mod(target_size(2),1)==0
                tmp=padarray(Pic_out1,[0 0 floor(size_diff(3)/2)],padval,'post');
                Pic_out=padarray(tmp,[0 0 ceil(size_diff(3)/2)],padval,'pre');
            else
                tmp=padarray(Pic_out1,[0 0 ceil(size_diff(3)/2)],padval,'post');
                Pic_out=padarray(tmp,[0 0 floor(size_diff(3)/2)],padval,'pre');
            end
        else
            if mod(size(Pic_in,3),2)==0
                vec=(ceil(abs(size_diff(3))/2)+1):ceil(abs(size_diff(3))/2)+target_size(3);
            else
                vec=(floor(abs(size_diff(3))/2)+1):floor(abs(size_diff(3))/2)+target_size(3);
            end
            Pic_out=Pic_out1(:,:,vec);
        end

    else %2D

        %pads images to a defined size; the frame around is filled with padval
        %targetsize can be a vector: [rows colmns]

        if length(target_size)==1
            target_size=[target_size target_size];
        end

        siz=size(Pic_in);
        size_diff=target_size(1:2)-siz(1:2);

        %adapt number of rows:
        if size_diff(1)>0 %if height of Pic_in is smaller than the target_size

            if mod(target_size(1),2)==0
                Pic_out1 = simplePadAboveBelow(Pic_in, [ceil(size_diff(1)/2), floor(size_diff(1)/2)], padval);
            else
                Pic_out1 = simplePadAboveBelow(Pic_in, [floor(size_diff(1)/2), ceil(size_diff(1)/2)], padval);
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
                Pic_out = simplePadLeftRight(Pic_out1, [ceil(size_diff(2)/2), floor(size_diff(2)/2)], padval);
            else
                Pic_out = simplePadLeftRight(Pic_out1, [floor(size_diff(2)/2), ceil(size_diff(2)/2)], padval);
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
end

function arrayOut = simplePadAboveBelow(arrayIn, nAboveBelow, padval)
    width = size(arrayIn, 2);
    arrayOut = [padval.*ones(nAboveBelow(1),width); arrayIn; padval.*ones(nAboveBelow(2),width)];
end

function arrayOut = simplePadLeftRight(arrayIn, nLeftRight, padval)
    height = size(arrayIn, 1);
    arrayOut = [padval.*ones(height, nLeftRight(1)), arrayIn, padval.*ones(height, nLeftRight(2))];
end