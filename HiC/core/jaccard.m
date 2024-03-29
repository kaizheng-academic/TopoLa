function jaccard_coefficient = jaccard(vector1, vector2)

intersection = length(intersect(vector1, vector2));
unionset = length(union(vector1, vector2));


jaccard_coefficient = intersection / unionset;
end
