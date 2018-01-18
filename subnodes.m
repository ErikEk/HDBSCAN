function children = subnodes(u,parents)

    % find children of this node
    subChildren = find(parents == u);
    if isempty( subChildren )
        children = [];
        return
    end
    
    children = subChildren;
    for child = subChildren
        newChildren = subnodes( child,parents );
        children = [children,newChildren];
    end
end