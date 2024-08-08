function [centroids, maxima] = findMaxima(image, slider)

    threshold = slider;
   
    maxima = imregionalmax(image);
    maxima = maxima & (image >= threshold);
    
    labeledMaxima = bwlabel(maxima);
    props = regionprops(labeledMaxima, 'Centroid');
    centroids = cat(1, props.Centroid);
    
    maxValues = [];
    if ~isempty(centroids)
        maxIndices = sub2ind(size(image), round(centroids(:,2)), round(centroids(:,1)));
        maxValues = image(maxIndices);
    end
       
    imshow(image, []);
    hold on
    if ~isempty(centroids)
        plot(centroids(:,1), centroids(:,2), 'r*');
    end
    hold off
end
