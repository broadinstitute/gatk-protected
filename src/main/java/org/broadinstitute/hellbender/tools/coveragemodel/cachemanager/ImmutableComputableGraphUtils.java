package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ImmutableComputableGraphUtils {

    /**
     * A simple builder class for {@link ImmutableComputableGraph}
     */
    public static class ImmutableComputableGraphBuilder {
        private final Set<CacheNode> nodes;
        private final Set<String> keys;
        private boolean cacheAutoUpdate;

        ImmutableComputableGraphBuilder() {
            nodes = new HashSet<>();
            keys = new HashSet<>();
            cacheAutoUpdate = false;
        }

        public ImmutableComputableGraphBuilder addPrimitiveNode(@Nonnull final String key,
                                                                @Nonnull final String[] tags,
                                                                @Nonnull Duplicable value) {
            Utils.nonNull(key);
            Utils.nonNull(tags);
            Utils.nonNull(value);
            Utils.validateArg(!keys.contains(key), "A node with key " + quote(key) + " already exists");
            nodes.add(new PrimitiveCacheNode(key, Arrays.stream(tags).collect(Collectors.toList()), value));
            keys.add(key);
            return this;
        }


        public ImmutableComputableGraphBuilder addNDArrayPrimitiveNode(@Nonnull final String key) {
            return addPrimitiveNode(key, new String[]{}, new DuplicableNDArray());
        }

        public ImmutableComputableGraphBuilder addComputableNode(@Nonnull final String key,
                                                                 @Nonnull final String[] tags,
                                                                 @Nonnull final String[] parents,
                                                                 @Nullable final ComputableNodeFunction func,
                                                                 final boolean cacheEvals) {
            Utils.nonNull(key);
            Utils.nonNull(tags);
            Utils.nonNull(parents);
            Utils.validateArg(!keys.contains(key), "A node with key " + quote(key) + " already exists");
            nodes.add(new ComputableCacheNode(key,
                    Arrays.stream(tags).collect(Collectors.toList()),
                    Arrays.stream(parents).collect(Collectors.toList()),
                    func, cacheEvals));
            keys.add(key);
            return this;
        }

        public ImmutableComputableGraphBuilder addExternallyComputableNode(@Nonnull final String key) {
            return addComputableNode(key, new String[] {}, new String[] {}, null, true);
        }

        public ImmutableComputableGraphBuilder enableCacheAutoUpdate() {
            cacheAutoUpdate = true;
            return this;
        }

        public ImmutableComputableGraphBuilder disableCacheAutoUpdate() {
            cacheAutoUpdate = false;
            return this;
        }

        public ImmutableComputableGraph build() {
            if (nodes.size() == 0) {
                throw new IllegalStateException("Can not make an empty cache node collection");
            } else {
                return new ImmutableComputableGraph(nodes, cacheAutoUpdate);
            }
        }
    }

    static String quote(final String str) {
        return "\"" + str + "\"";
    }
}
