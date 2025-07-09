window.ega_curve_mo = function (el, data, opts) {
  const margin = { top: 50, right: 20, bottom: 70, left: 50 };
  const chartContainer = document.querySelector(el).closest(".chart-container");
  const svgContainer = document.querySelector(el).closest(".svg-container");
  const svg = d3.select(el);
  let isZoomActive = false;

  function debounce(cb, delay = 1000) {
    let timeout;
    return (...args) => {
      clearTimeout(timeout);
      timeout = setTimeout(() => cb(...args), delay);
    };
  }
  function getActiveChartValue(container){
    const checkedArray = container.querySelector("input[type='radio']:checked");
    return checkedArray.value;
  }
  function getSVGDimensions(){
    let width = opts.width || 680;
    let height = opts.height || 253;
    if(opts.width ==="auto"){       
      const activeChart = chartContainer.querySelector(
        `section.wide[value="${getActiveChartValue(chartContainer)}"]`
      );
        
      const svgContainerDimensions = activeChart.getBoundingClientRect();
      width = svgContainerDimensions.width;
      height = 488;
    }
    return [width,height]
  }
  function getSeparators(data){    
    const arr = [];
    let current = null;
    let labelsArray = [];
    for (const dataPoint of data) {
      if(!current || !labelsArray.includes(dataPoint[0])){
        current = dataPoint;
        labelsArray.push(dataPoint[0])
        arr.push(dataPoint);
      }
    }    
    return arr;
  }
  function expandScaleDomain(scale,amount){
    const newDomain = scale.domain();
    if(typeof amount === 'string' && amount.endsWith("px")){
      const amount_px = Number(amount.replace("px",""));
      const amount_domain = scale.invert(amount_px);
      newDomain[1]+=amount_domain;
    }else{
      const ticks = scale.ticks();
      const scaleStep = ticks[1] - ticks[0];
      newDomain[1]+=amount*scaleStep;
    }
    scale.domain(newDomain);
  }
  function addTooltip(x,y,svg,data,opts){
    const tooltip = svg.append("g");
    // Add the event listeners that show or hide the tooltip.
    const formatX = d3.format('~s')
    const bisect = d3.bisector(d => d[1]).center;

    svg
      .on("pointerenter pointermove", pointermoved)
      .on("pointerleave", pointerleft)
      .on("touchstart", event => event.preventDefault());


    function pointermoved(event) {
      const i = bisect(data, x.invert(d3.pointer(event)[0]));
      const cx = x(data[i][1]);
      const cy = y(data[i][2]);
      tooltip.style("display", null);
      tooltip.attr("transform", `translate(${cx},${cy})`);

      const path = tooltip.selectAll("path")
        .data([,])
        .join("path")
          .attr("stroke", "black")
          .attr("class","graph-tooltip")

      const text = tooltip.selectAll("text")
        .data([,])
        .join("text")
        .call(text => text
          .selectAll("tspan")
          .data([opts.title.toLowerCase(),`${formatX(data[i][1])}: ${data[i][2]}`])
          .join("tspan")
            .attr("x", 0)
            .attr("y", (_, i) => `${i * 1.1}em`)
            .attr("font-weight", (_, i) => i ? null : "bold")
            .text(d => d));

      // Measure text and decide if it needs to flip
      const { width: w } = text.node().getBBox();
      const padding = 20;
      const flipped = (cx + w + padding + 20) > svg.node().getBoundingClientRect().width;    

      size(text, path, flipped);
    }

    function pointerleft() {
      tooltip.style("display", "none");
    }
    
    // Wraps the text with a callout path of the correct size, as measured in the page.
    function size(text, path, flipped = false) {
      const { x, y, width: w, height: h } = text.node().getBBox();
      const paddingX = 15;
      const offsetY = 12;
      if (flipped) {
        text.attr("transform", `translate(${-w - paddingX},${offsetY - h / 2})`);
        path.attr("d", `M${-5},${-h / 2 - 3} V-5 L0,0 L-5,5 V${h / 2 + 3} H${-w - 2 * paddingX} V${-h / 2 - 3} Z`);
      } else {
        text.attr("transform", `translate(${paddingX - x},${offsetY - h / 2})`);
        path.attr("d", `M5,${-h / 2 - 3} V-5 L0,0 L5,5 V${h / 2 + 3} H${w + 2 * paddingX} V${-h / 2 - 3} Z`);
      }
    }
  }
  function getResizeObserver(){
    let chartSize = 0;
    let loading = false;
    const svgNode = d3.select(el).node();

    const debouncedBuildChart = debounce((newSize) => {
    chartSize = newSize;
    svgContainer.classList.remove("loading");
    loading = false;
    // console.log("building chart...");
    buildChart();
    }, 50);

    const observer = new ResizeObserver(([entry]) => {
    if (!entry.borderBoxSize) return;

    const newSize = entry.borderBoxSize[0].inlineSize;
    if (newSize === chartSize) return;

    if (!loading) {
      svgContainer.classList.add("loading");
      svgNode.innerHTML = "";
      loading = true;
      // console.log("adding loader...");
    }

    // Update chart size and rebuild after debounce delay
    debouncedBuildChart(newSize);
    });
    return observer;
  }
  function updateChart(event, x) {
    const dx = event.selection;
    if (dx) {
      const [from, to] = [
        x.invert(dx[0]) > 0 ? x.invert(dx[0] < 60 && x.domain()[0] === 0 ? 0 : dx[0]) : 0,
        x.invert(dx[1]),
      ];        
      d3.select(el).node().innerHTML="";
      buildChart(from,to);
      isZoomActive = true;
    }
  }
  function resetZoom(){
    if (isZoomActive) {
       isZoomActive = false;
       d3.select(el).node().innerHTML="";
       buildChart();
     }
  }
  function buildChart(from = 0, to = Infinity){
    const [width,height] = getSVGDimensions();
    const filteredData =
      from === 0 && to === Infinity
        ? data
        : data.filter((point) => point[1] >= from && point[1] < to);
          
    // Declare the x (horizontal position) scale.
    const x = d3.scaleLinear(d3.extent(filteredData, d => d[1]), [margin.left, width - margin.right]);
    // Declare the y (vertical position) scale.
    const y = d3.scaleLinear([0, d3.max(filteredData, d => d[2])], [height - margin.bottom, margin.top]).nice();

    //expand scale domain (making extra space in the graph)
    if(from === 0 && to === Infinity){
      if(opts.expandX) expandScaleDomain(x,opts.expandX)
      if(opts.expandY) expandScaleDomain(y,opts.expandY)
    }else{
      if(opts.expandXOnZoom) expandScaleDomain(x,opts.expandXOnZoom)
      if(opts.expandYOnZoom) expandScaleDomain(y,opts.expandYOnZoom)
    }
    // Declare the line generator.
    const line = d3.line()
    .x(d => x(d[1]))
    .y(d => y(d[2]));

    // Declare separators
    const separators = getSeparators(filteredData);
    

    svg
     .attr("width", width)
     .attr("height", height)
     .attr("viewBox", [0, 0, width, height])
     .attr("style", "max-width: 100%; height: auto; height: intrinsic;")

    // Add the x-axis.
    svg.append("g")
    .attr("class", "axis x-axis")
    .attr("transform", `translate(0,${height - margin.bottom})`)
    .call(d3.axisBottom(x).ticks(10,d3.format('~s')).tickSizeOuter(0))
    .call(g => g.selectAll(".tick line").clone()
    .attr("y2", -(height - margin.top - margin.bottom))
    .attr("stroke-opacity", 0.1));


    // Add the y-axis, add grid lines.
    svg.append("g")
    .attr("class", "axis y-axis")
    .attr("transform", `translate(${margin.left},0)`)
    .call(d3.axisLeft(y).ticks(10, d3.format('~s')))
    .call(g => g.selectAll(".tick line").clone()
        .attr("x2", width - margin.left - margin.right)
        .attr("stroke-opacity", 0.1));   


    // Append a path for the line.
    svg.append("path")
    .attr("class", "plot")
    .attr('d', line(filteredData));

    // x-axis name
    if(opts.xAxis !== undefined){
      svg.append("text")
          .attr("class", 'x-axis-name')
          .attr("x", width/2)
          .attr("y", height - (margin.bottom/4))
          .attr("text-anchor", "middle")
          .text(opts.xAxis);
    }
    // y-axis name
    if(opts.yAxis !== undefined ){
      svg.append("text")
          .attr("class", 'y-axis-name')
          .attr("x", -height/2 -10)
          .attr("y", margin.left/4 - 7)
          .attr("dy", "0.35em")
          .attr("text-anchor", "middle")
          .attr("transform", "rotate(-90)")
          .text(opts.yAxis);
    }
    //title
    if(opts.showTitle && opts.title !== undefined){
      svg.append("text")
          .attr("class", 'title')
          .attr("x", width/2)
          .attr("y", (margin.top/3))
          .attr("text-anchor", "middle")
          .text(opts.title);
    }

    const separatorGroup = svg.append("g").attr("class", "separators");

    // Append separators
    separatorGroup.selectAll(".separator-line") //selects the existent separators
    .data(separators.filter((d,i)=>i>0)) //binds them to data
    .enter() //selects the other needed elements (the one data requires but are not in the previous selection)
    .append("line") //creates a line for each element
    .attr("class", "separator-line")
    .attr("x1", d => x(d[1])) // x position
    .attr("x2", d => x(d[1]))
    .attr("y1", margin.top) // Start of the vertical line
    .attr("y2", height-margin.bottom) // End of the vertical line
    .attr("stroke", "black")
    .attr("stroke-dasharray", "4,4"); // Optional dashed line

    // Append labels
    separatorGroup.selectAll(".separator-label")
    .data(separators)
    .enter()
    .append("text")
    .attr("class", "separator-label")
    .attr("transform", "rotate(-90)")
    .attr("y", d => x(d[1]) + 7) // x pos
    .attr("x", -1*(margin.top + 35)) // y pos
    .attr("dy", "0.35em")
    .text(d => d[0])
    .attr("font-size", "12px")
    .attr("fill", "black");

    // Add tooltip
    if(opts.showTooltip)
      addTooltip(x,y,svg,filteredData,opts);
    
    // Add brush zoom
    const brush = d3
      .brushX()
      .extent([
        [margin.left, margin.top],
        [width - margin.right , height - margin.bottom],
      ])
      .on("end", (event) => updateChart(event,x));

    svg.append("g").attr("class", "brush").call(brush);
    svg.on("dblclick", resetZoom);
  }

  const resizeObserver = getResizeObserver();
  resizeObserver.observe(chartContainer);

};