

function model=CreateModel()

    % Source
    xs=0;
    ys=0;
    
    % Target (Destination)
    xt=9;
    yt=8;
    
    xobs=[1.5 6.0 1.2 8.1];
    yobs=[4.5 7.0 1.5 1.2];
    robs=[1.5 0.7 0.8 1.2];
    
    n=10000;
    
    xmin=-10;
    xmax= 10;
    
    ymin=-10;
    ymax= 10;
    
    model.xs=xs;
    model.ys=ys;
    model.xt=xt;
    model.yt=yt;
    model.xobs=xobs;
    model.yobs=yobs;
    model.robs=robs;
    model.n=n;
    model.xmin=xmin;
    model.xmax=xmax;
    model.ymin=ymin;
    model.ymax=ymax;
    
end