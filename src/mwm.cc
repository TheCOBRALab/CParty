#include "mwm.hh"

#include <queue>
#include <cassert>
#include <iostream>

maxWeightMatching::maxWeightMatching(std::vector<elem_prob_s> &edges,cand_pos_t n) : n(n), nedge(edges.size()), mate(n,-1), label(2*n,0),endpoint(2*nedge),neighbend(n),labelend(2*n,-1),inblossom(n),blossomparent(2*n,-1),blossomchilds(2*n),blossombase(2*n,-1),blossomendps(2*n),bestedge(2*n,-1),blossombestedges(2*n),unusedblossoms(n),dualvar(2*n,0),allowedge(nedge,0){
    this->edges = edges;
    // Find the maximum edge weight.
    maxweight = std::max(0.,std::max_element(edges.begin(),edges.end(),[](elem_prob_s a, elem_prob_s b) {return a.p < b.p;})->p);
    // Probably a faster way without division
    for(cand_pos_t i = 0; i<2*nedge;++i){
        if(i%2 == 0){
            endpoint[i] = edges[i/2].i;
        }
        else{
            endpoint[i] = edges[i/2].j;
        }
    }
    for (cand_pos_t k =0;k<nedge;++k){
        elem_prob_s *e = &edges[k];
        neighbend[e->i].push_back(2*k+1);
        neighbend[e->j].push_back(2*k);
    }
    for(cand_pos_t i = 0; i<n;++i) inblossom[i] = i;
    for(cand_pos_t i = 0; i<n;++i) blossombase[i] = i;
    for(cand_pos_t i = n; i<2*n;++i) unusedblossoms[i-n] = i; //better way to do this?
    for(cand_pos_t i = 0; i<n;++i) dualvar[i] = maxweight;
    
}


pf_t maxWeightMatching::slack(cand_pos_t k){
    elem_prob_s *e = &edges[k];
    return dualvar[e->i] + dualvar[e->j] - 2 * e->p;
}

void maxWeightMatching::blossomLeaves(cand_pos_t b,std::vector<cand_pos_t> &leaves){
    if(b < n){
        leaves.push_back(b);
    } else{
        for(cand_pos_t t:blossomchilds[b]){
            if(t < n){
                leaves.push_back(t);
            } else{
                blossomLeaves(t,leaves);
            }
        }
    }
}

void maxWeightMatching::assignLabel(cand_pos_t w, cand_pos_t t, cand_pos_t p){
    cand_pos_t b = inblossom[w];
    assert(label[w] == 0 && label[b] == 0);
    label[w] = label[b] = t;
    labelend[w] = labelend[b] = p;
    bestedge[w] = bestedge[b] = -1;
    if(t == 1){
        std::vector<cand_pos_t> leaves;
        blossomLeaves(b, leaves);
        // b became an S-vertex/blossom; add it(s vertices) to the queue.
        for(cand_pos_t v: leaves) queue.push_back(v);
    }
    else if(t == 2){
        // b became a T-vertex/blossom; assign label S to its mate.
        // (If b is a non-trivial blossom, its base is the only vertex
        // with an external mate.)
        cand_pos_t base = blossombase[b];
        assert(mate[base] >= 0);
        assignLabel(endpoint[mate[base]], 1, mate[base] ^ 1);
    }
}

cand_pos_t maxWeightMatching::scanBlossom(cand_pos_t v, cand_pos_t w){ //This maybe should be references, maybe not
    // Trace back from v and w, placing breadcrumbs as we go.
    std::vector<cand_pos_t> path;
    cand_pos_t base = -1;
    cand_pos_t b;
    while(v != -1 or w != -1){
        // Look for a breadcrumb in v's blossom or put a new breadcrumb.
        b = inblossom[v];
        if(label[b] & 4){
            base = blossombase[b];
            break;
        }
        assert(label[b] == 1);
        path.push_back(b);
        label[b] = 5;
        // Trace one step back.
        assert(labelend[b] == mate[blossombase[b]]);
        if(labelend[b] == -1){
            // The base of blossom b is single; stop tracing this path.
            v = -1;
        } else{
            v = endpoint[labelend[b]];
            b = inblossom[v];
            assert(label[b] == 2);
            // b is a T-blossom; trace one more step back.
            assert(labelend[b] >= 0);
            v = endpoint[labelend[b]];
        }
        // Swap v and w so that we alternate between both paths.
        if(w != -1) std::swap(v,w);
    }
    // Remove breadcrumbs.
    for(cand_pos_t b : path) label[b] = 1;
    // Return base vertex, if we found one.
    return base;
}

void maxWeightMatching::addBlossom(cand_pos_t base, cand_pos_t k){
    elem_prob_s *e = &edges[k]; //i is v, j is w, p is wt
    cand_pos_t v = e->i;
    cand_pos_t w = e->j;
    cand_pos_t bb = inblossom[base];
    cand_pos_t bv = inblossom[v];
    cand_pos_t bw = inblossom[w];
    // Create blossom.
    cand_pos_t b = unusedblossoms.back();
    unusedblossoms.pop_back();
    blossombase[b] = base;
    blossomparent[b] = -1;
    blossomparent[bb] = b;
    // Make list of sub-blossoms and their interconnecting edge endpoints.
    blossomchilds[b].clear();
    std::vector<cand_pos_t> path;
    blossomendps[b].clear();
    std::vector<cand_pos_t> endps;
    // Trace back from v to base.
    while(bv != bb){
        // Add bv to the new blossom.
        blossomparent[bv] = b;
        path.push_back(bv);
        endps.push_back(labelend[bv]);
        assert(label[bv] == 2 || (label[bv] == 1 && labelend[bv] == mate[blossombase[bv]]));
        // Trace one step back.
        assert(labelend[bv] >= 0);
        v = endpoint[labelend[bv]];
        bv = inblossom[v];
    }
    // Reverse lists, add endpoint that connects the pair of S vertices.
    path.push_back(bb);
    reverse(path.begin(), path.end());
    reverse(endps.begin(), endps.end());
    endps.push_back(2*k);
    // Trace back from w to base.
    while(bw != bb){
        // Add bw to the new blossom.
        blossomparent[bw] = b;
        path.push_back(bw);
        endps.push_back(labelend[bw] ^ 1);
        assert (label[bw] == 2 || (label[bw] == 1 && labelend[bw] == mate[blossombase[bw]]));
        // Trace one step back.
        assert(labelend[bw] >= 0);
        w = endpoint[labelend[bw]];
        bw = inblossom[w];
    }
    // Set label to S.
    assert(label[bb] == 1);
    label[b] = 1;
    labelend[b] = labelend[bb];
    // Set dual variable to zero.
    dualvar[b] = 0;
    // Relabel vertices.
    std::vector<cand_pos_t> leavesb;
    blossomLeaves(b,leavesb);
    for(cand_pos_t v : leavesb){
        if(label[inblossom[v]] == 2)
            // This T-vertex now turns into an S-vertex because it becomes
            // part of an S-blossom; add it to the queue.
            queue.push_back(v);
        inblossom[v] = b;
    }
    // Compute blossombestedges[b].
    std::vector<cand_pos_t> bestedgeto(2*n,-1);
    for(cand_pos_t bv : path){
        std::vector<std::vector<cand_pos_t> > nblists;
        if(blossombestedges[bv].empty()){
            // This subblossom does not have a list of least-slack edges;
            // get the information from the vertices.
            std::vector<cand_pos_t> leavesv;
            blossomLeaves(bv,leavesv);
            for(cand_pos_t i = 0;i<(cand_pos_t) leavesv.size();++i){
                std::vector<cand_pos_t> temp;
                for(cand_pos_t j = 0;j<(cand_pos_t) neighbend[leavesv[i]].size();++j){
                    temp.push_back(neighbend[leavesv[i]][j]);
                }
                nblists.push_back(temp);
            }
        } else{
            // Walk this subblossom's least-slack edges.
            nblists.push_back(blossombestedges[bv]);
        }
        for(std::vector<cand_pos_t> nblist:nblists){
            for(cand_pos_t k:nblist){
                cand_pos_t i = edges[k].i;
                cand_pos_t j = edges[k].j;
                if(inblossom[j] == b) std::swap(i,j);
                cand_pos_t bj = inblossom[j];
                if (bj != b && label[bj] == 1 && (bestedgeto[bj] == -1 || slack(k) < slack(bestedgeto[bj])))
                    bestedgeto[bj] = k;
            }
        }
        // Forget about least-slack edges of the subblossom.
        blossombestedges[bv].clear();
        bestedge[bv] = -1;
    }
    blossombestedges[b].clear();
    for (cand_pos_t k : bestedgeto) {
        if (k != -1) blossombestedges[b].push_back(k);
    }
    // Select bestedge[b].
    bestedge[b] = -1;
    for(cand_pos_t k : blossombestedges[b])
        if(bestedge[b] == -1 || slack(k) < slack(bestedge[b])) bestedge[b] = k;
}

void maxWeightMatching::expandBlossom(cand_pos_t b, bool endstage){
    // Convert sub-blossoms into top-level blossoms.
    for(cand_pos_t s:blossomchilds[b]){
        blossomparent[s] = -1;
        if(s < n){
            inblossom[s] = s;
        } else if(endstage && dualvar[s] == 0){
            // Recursively expand this sub-blossom.
            expandBlossom(s, endstage);
        } else{
            std::vector<cand_pos_t> leavesv;
            blossomLeaves(s,leavesv);
            for(cand_pos_t v : leavesv)
                inblossom[v] = s;
        }
    }
    // If we expand a T-blossom during a stage, its sub-blossoms must be
    // relabeled.
    if (!endstage && label[b] == 2){
        // Start at the sub-blossom through which the expanding
        // blossom obtained its label, and relabel sub-blossoms untili
        // we reach the base.
        // Figure out through which sub-blossom the expanding blossom
        // obtained its label initially.
        assert(labelend[b] >= 0);
        cand_pos_t entrychild = inblossom[endpoint[labelend[b] ^ 1]];
        // Decide in which direction we will go round the blossom.
        auto it = std::find(blossomchilds[b].begin(), blossomchilds[b].end(), entrychild);
        cand_pos_t j = std::distance(blossomchilds[b].begin(), it);
        cand_pos_t jstep,endptrick;
        if(j & 1){
            // Start index is odd; go forward and wrap.
            j -= blossomchilds[b].size();
            jstep = 1;
            endptrick = 0;
        }
        else{
            // Start index is even; go backward.
            jstep = -1;
            endptrick = 1;
        }
        // Move along the blossom until we get to the base.
        cand_pos_t p = labelend[b];
        while(j != 0){
            // Relabel the T-sub-blossom.
            label[endpoint[p ^ 1]] = 0;
            label[endpoint[blossomendps[b][j-endptrick]^endptrick^1]] = 0;
            assignLabel(endpoint[p ^ 1], 2, p);
            // Step to the next S-sub-blossom and note its forward endpoint.
            allowedge[blossomendps[b][j-endptrick]/2] = true;
            j += jstep;
            p = blossomendps[b][j-endptrick] ^ endptrick;
            // Step to the next T-sub-blossom.
            allowedge[p/2] = true;
            j += jstep;
        }
        // Relabel the base T-sub-blossom WITHOUT stepping through to
        // its mate (so don't call assignLabel).
        cand_pos_t bv = blossomchilds[b][j];
        label[endpoint[p ^ 1]] = label[bv] = 2;
        labelend[endpoint[p ^ 1]] = labelend[bv] = p;
        bestedge[bv] = -1;
        // Continue along the blossom until we get back to entrychild.
        j += jstep;
        while(blossomchilds[b][j] != entrychild){
            // Examine the vertices of the sub-blossom to see whether
            // it is reachable from a neighbouring S-vertex outside the
            // expanding blossom.
            cand_pos_t bv = blossomchilds[b][j];
            if(label[bv] == 1){
                // This sub-blossom just got label S through one of its
                // neighbours; leave it.
                j += jstep;
                continue;
            }
            cand_pos_t v = 0;
            std::vector<cand_pos_t> leavesbv;
            blossomLeaves(bv,leavesbv);
            for(cand_pos_t bvv: leavesbv)
                if(label[bvv] != 0) {v=bvv; break;}
            // If the sub-blossom contains a reachable vertex, assign
            // label T to the sub-blossom.
            if(label[v] != 0){
                assert(label[v] == 2);
                assert(inblossom[v] == bv);
                label[v] = 0;
                label[endpoint[mate[blossombase[bv]]]] = 0;
                assignLabel(v, 2, labelend[v]);
            }
            j += jstep;
        }
    }
    // Recycle the blossom number.
    label[b] = labelend[b] = -1;
    blossomchilds[b].clear(); blossomendps[b].clear();
    blossombase[b] = -1;
    blossombestedges[b].clear();
    bestedge[b] = -1;
    unusedblossoms.push_back(b);
}

void maxWeightMatching::augmentBlossom(cand_pos_t b, cand_pos_t v){
    // Bubble up through the blossom tree from vertex v to an immediate
    // sub-blossom of b.
    cand_pos_t t = v;
    cand_pos_t jstep,endptrick;
    while(blossomparent[t] != b) t = blossomparent[t];

    // Recursively deal with the first sub-blossom.
    if (t >= n) augmentBlossom(t, v);
    // Decide in which direction we will go round the blossom.
    cand_pos_t i,j;
    auto it = std::find(blossomchilds[b].begin(), blossomchilds[b].end(), t);
    i = j = std::distance(blossomchilds[b].begin(), it);
    if(i & 1){
        // Start index is odd; go forward and wrap.
        j -= blossomchilds[b].size();
        jstep = 1;
        endptrick = 0;
    } else{
        // Start index is even; go backward.
        jstep = -1;
        endptrick = 1;
    }
    // Move along the blossom until we get to the base.
    while(j != 0){
        // Step to the next sub-blossom and augment it recursively.
        j += jstep;;
        t = blossomchilds[b][j];
        cand_pos_t p = blossomendps[b][j-endptrick] ^ endptrick;
        if(t >= n) augmentBlossom(t, endpoint[p]);
        // Step to the next sub-blossom and augment it recursively.
        j += jstep;
        t = blossomchilds[b][j];
        if(t >= n) augmentBlossom(t, endpoint[p ^ 1]);
        // Match the edge connecting those sub-blossoms.
        mate[endpoint[p]] = p ^ 1;
        mate[endpoint[p ^ 1]] = p;
    }
    // Rotate the list of sub-blossoms to put the new base at the front.
    std::rotate(blossomchilds[b].begin(),blossomchilds[b].begin() + i,blossomchilds[b].end());
    std::rotate(blossomendps[b].begin(),blossomendps[b].begin() + i,blossomendps[b].end());
    blossombase[b] = blossombase[blossomchilds[b][0]];
    assert(blossombase[b] == v);
}

void maxWeightMatching::augmentMatching(cand_pos_t k){
    cand_pos_t v = edges[k].i;
    cand_pos_t w = edges[k].j;
    // for (s, p) in ((v, 2*k+1), (w, 2*k)){
    for (cand_pos_t pass = 0; pass < 2; ++pass) {
        cand_pos_t s = (pass == 0 ? v : w);
        cand_pos_t p = (pass == 0 ? 2*k+1 : 2*k);
        // Match vertex s to remote endpoint p. Then trace back from s
        // until we find a single vertex, swapping matched and unmatched
        // edges as we go.
        while(true){
            cand_pos_t bs = inblossom[s];
            assert(label[bs] == 1);
            assert(labelend[bs] == mate[blossombase[bs]]);
            // Augment through the S-blossom from s to base.
            if(bs >= n) augmentBlossom(bs, s);
            // Update mate[s]
            mate[s] = p;
            // Trace one step back.
            if(labelend[bs] == -1)
                // Reached single vertex; stop.
                break;
            cand_pos_t t = endpoint[labelend[bs]];
            cand_pos_t bt = inblossom[t];
            assert(label[bt] == 2);
            // Trace one step back.
            assert(labelend[bt] >= 0);
            s = endpoint[labelend[bt]];
            cand_pos_t j = endpoint[labelend[bt] ^ 1];
            // Augment through the T-blossom from j to base.
            assert(blossombase[bt] == t);
            if(bt >= n) augmentBlossom(bt, j);
            // Update mate[j]
            mate[j] = labelend[bt];
            // Keep the opposite endpoint;
            // it will be assigned to mate[s] in the next step.
            p = labelend[bt] ^ 1;
        }
    }
}

pf_t maxWeightMatching::maximumWeightedMatching(std::vector<pf_t> &pu,std::string &structure, sparse_tree &tree){


    // Main loop: continue until no further improvement is possible.
    for(cand_pos_t t=0;t<n;++t){

        // Each iteration of this loop is a "stage".
        // A stage finds an augmenting path and uses that to improve
        // the matching.

        // Remove labels from top-level blossoms/vertices.
        label = std::vector<cand_pos_t>(2*n,0);

        // Forget all about least-slack edges.
        bestedge = std::vector<cand_pos_t>(2*n,-1);
        for(cand_pos_t i = n;i<2*n;++i) blossombestedges.clear();

        // Loss of labeling means that we can not be sure that currently
        // allowable edges remain allowable througout this stage.
        allowedge = std::vector<char>(nedge,false);

        // Make queue empty.
        queue.clear();

        // Label single blossoms/vertices with S and put them in the queue.
        for(cand_pos_t v=0;v<n;++v) if(mate[v] == -1 && label[inblossom[v]] == 0) assignLabel(v, 1, -1);

        // Loop until we succeed in augmenting the matching.
        bool augmented = 0;
        while(true){

            // Each iteration of this loop is a "substage".
            // A substage tries to find an augmenting path;
            // if found, the path is used to improve the matching and
            // the stage ends. If there is no augmenting path, the
            // primal-dual method is used to pump some slack out of
            // the dual variables.

            // Continue labeling until all vertices which are reachable
            // through an alternating path have got a label.
            while(!queue.empty() && !augmented){

                // Take an S vertex from the queue.
                cand_pos_t v = queue.back();
                queue.pop_back();
                assert(label[inblossom[v]] == 1);

                // Scan its neighbours:
                for(cand_pos_t p : neighbend[v]){
                    cand_pos_t k = p / 2;
                    cand_pos_t w = endpoint[p];
                    pf_t kslack = 0;
                    // w is a neighbour to v
                    if(inblossom[v] == inblossom[w])
                        // this edge is internal to a blossom; ignore it
                        continue;
                    if(!allowedge[k]){
                        kslack = slack(k);
                        if(kslack <= 0){
                            // edge k has zero slack => it is allowable
                            allowedge[k] = true;
                        }
                    }
                    if(allowedge[k]){
                        if(label[inblossom[w]] == 0){
                            // (C1) w is a free vertex;
                            // label w with T and label its mate with S (R12).
                            assignLabel(w, 2, p ^ 1);
                        } else if(label[inblossom[w]] == 1){
                            // (C2) w is an S-vertex (not in the same blossom);
                            // follow back-links to discover either an
                            // augmenting path or a new blossom.
                            cand_pos_t base = scanBlossom(v, w);
                            if(base >= 0){
                                // Found a new blossom; add it to the blossom
                                // bookkeeping and turn it into an S-blossom.
                                addBlossom(base, k);
                            } else{
                                // Found an augmenting path; augment the
                                // matching and end this stage.
                                augmentMatching(k);
                                augmented = 1;
                                break;
                            }
                        } else if(label[w] == 0){
                            // w is inside a T-blossom, but w itself has not
                            // yet been reached from outside the blossom;
                            // mark it as reached (we need this to relabel
                            // during T-blossom expansion).
                            assert(label[inblossom[w]] == 2);
                            label[w] = 2;
                            labelend[w] = p ^ 1;
                        }
                    } else if(label[inblossom[w]] == 1){
                        // keep track of the least-slack non-allowable edge to
                        // a different S-blossom.
                        cand_pos_t b = inblossom[v];
                        if(bestedge[b] == -1 or kslack < slack(bestedge[b])){
                            bestedge[b] = k;
                        }
                    } else if(label[w] == 0){
                        // w is a free vertex (or an unreached vertex inside
                        // a T-blossom) but we can not reach it yet;
                        // keep track of the least-slack edge that reaches w.
                        if(bestedge[w] == -1 || kslack < slack(bestedge[w]))
                            bestedge[w] = k;
                    }
                }
            }
            if(augmented) break;

            // There is no augmenting path under these constraints;
            // compute delta and reduce slack in the optimization problem.
            // (Note that our vertex dual variables, edge slacks and delta's
            // are pre-multiplied by two.)
            cand_pos_t deltatype = -1;
            cand_pos_t delta,deltaedge,deltablossom;

            // Compute delta1: the minumum value of any vertex dual.
            deltatype = 1;
            delta = *std::min_element(dualvar.begin(), dualvar.begin()+n); //min(dualvar[:n]);

            // Compute delta2: the minimum slack on any edge between
            // an S-vertex and a free vertex.
            for(cand_pos_t v=0;v<n;++v){
                if(label[inblossom[v]] == 0 and bestedge[v] != -1){
                    cand_pos_t d = slack(bestedge[v]);
                    if(deltatype == -1 or d < delta){
                        delta = d;
                        deltatype = 2;
                        deltaedge = bestedge[v];
                    }
                }
            }

            // Compute delta3: half the minimum slack on any edge between
            // a pair of S-blossoms.
            for(cand_pos_t b=0;b<2*n;++b){
                if (blossomparent[b] == -1 && label[b] == 1 && bestedge[b] != -1 ){
                    pf_t kslack = slack(bestedge[b]);
                    cand_pos_t d = kslack / 2;
                    if(deltatype == -1 or d < delta){
                        delta = d;
                        deltatype = 3;
                        deltaedge = bestedge[b];
                    }
                }
            }

            // Compute delta4: minimum z variable of any T-blossom.
            for(cand_pos_t b=n;b<2*n;++b){
                if ( blossombase[b] >= 0 && blossomparent[b] == -1 && label[b] == 2 && (deltatype == -1 || dualvar[b] < delta) ){
                    delta = dualvar[b];
                    deltatype = 4;
                    deltablossom = b;
                }
            }

            // Update dual variables according to delta.
            for(cand_pos_t v=0;v<n;++v){
                if(label[inblossom[v]] == 1){
                    // S-vertex: 2*u = 2*u - 2*delta
                    dualvar[v] -= delta;
                } else if(label[inblossom[v]] == 2){
                    // T-vertex: 2*u = 2*u + 2*delta
                    dualvar[v] += delta;
                }
            }
            for(cand_pos_t b=n;b<2*n;++b){
                if(blossombase[b] >= 0 && blossomparent[b] == -1){
                    if(label[b] == 1){
                        // top-level S-blossom: z = z + 2*delta
                        dualvar[b] += delta;
                    } else if(label[b] == 2){
                        // top-level T-blossom: z = z - 2*delta
                        dualvar[b] -= delta;
                    }
                }
            }

            // Take action at the point where minimum delta occurred.
            if(deltatype == 1){ 
                // No further improvement possible; optimum reached.
                break;
            } else if(deltatype == 2){
                // Use the least-slack edge to continue the search.
                allowedge[deltaedge] = true;
                cand_pos_t i = edges[deltaedge].i;
                cand_pos_t j = edges[deltaedge].j;
                if(label[inblossom[i]] == 0) std::swap(i,j);
                assert(label[inblossom[i]] == 1);
                queue.push_back(i);
            } else if(deltatype == 3){
                // Use the least-slack edge to continue the search.
                allowedge[deltaedge] = true;
                cand_pos_t i = edges[deltaedge].i;
                assert(label[inblossom[i]] == 1);
                queue.push_back(i);
            }
            else if(deltatype == 4){
                // Expand the least-z blossom.
                expandBlossom(deltablossom, false);
            }
        }
        // Stop when no more augmenting path can be found.
        if(!augmented) break;

        // End of a stage; expand all S-blossoms which have dualvar = 0.
        for (cand_pos_t b =n;b<2*n;++b)
            if ( blossomparent[b] == -1 && blossombase[b] >= 0 && label[b] == 1 && dualvar[b] == 0 ) expandBlossom(b, true);
    }
    // Transform mate[] such that mate[v] is the vertex to which v is paired.
    for(cand_pos_t v =0;v<n;++v) if(mate[v] >= 0) mate[v] = endpoint[mate[v]];
    for(cand_pos_t v =0;v<n;++v) assert(mate[v] == -1 or mate[mate[v]] == v);

    pf_t MEA=0.0;

    for(elem_prob_s &e: edges){
        if(mate[e.i]==e.j) {
            MEA+=2*e.p;
            if(tree.weakly_closed(e.i,e.j)){
                structure[e.i-1] = '(';
                structure[e.j-1] = ')';
            } else{
                structure[e.i-1] = '[';
                structure[e.j-1] = ']';
            }
        }
    }
    for(cand_pos_t i = n; i>0;--i) if(mate[i] == -1) MEA += pu[i];

    // return mate;
    return MEA;
}
